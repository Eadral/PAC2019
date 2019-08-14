module solver_pcg_mod

   use kinds_mod
   use netcdf_mod
   use blocks
   use distribution
   use domain
   use domain_size
   use constants
   use boundary
   use global_reductions
   use gather_scatter
   use broadcast
   use grid
   use io
   use time_management
   use exit_mod
   use communicate, only: my_task, master_task
   implicit none
   save
   public :: pcg
   real (r8), dimension (nx_block,ny_block,max_blocks_clinic) :: &
      AC,                &! time-independent part of center 9pt weight
      A0_CLINIC           ! time-dependent center weight of 9pt operator
                          !   in baroclinic block distribution

   integer (int_kind) :: &
      solv_sum_iters      ! accumulated no of iterations (diagnostic)

   real (r8) ::  &
      rms_residual        ! residual (also a diagnostic)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  other operator and preconditioner weights for barotropic operator
!
!-----------------------------------------------------------------------

   real (r8), dimension (nx_block,ny_block,max_blocks_tropic) :: &
      A0,AN,AE,ANE,         &! barotropic (9pt) operator coefficients
      RCALCT_B               ! land mask in barotropic distribution

   real (r8), dimension (:,:,:), allocatable :: &
      PCC,PCN,PCS,PCE,PCW,  &! preconditioner coefficients
      PCNE,PCSE,PCNW,PCSW

!-----------------------------------------------------------------------
!
!  scalar convergence-related variables
!
!-----------------------------------------------------------------------

   logical (log_kind) :: &
      lprecond            ! true if computed preconditioner to be used

   real (r8) ::          &
      solv_convrg,       &! convergence error criterion
      sor,               &! for jacobi solver
      resid_norm          ! residual normalization

   integer (int_kind), parameter :: &
      solv_pcg = 1,      &! predefined solver types
      solv_cgr = 2,      &
      solv_jac = 3

   integer (int_kind) :: &
      solv_itype,        &! integer solver method (1=pcg, 2=cgr, 3=jac)
      solv_max_iters,    &! max number of iterations for solver
      solv_ncheck         ! check convergence every ncheck iterations

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      B                         ! right hand side of linear system

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      X,X0                  ! on input,  an initial guess for the solution
                         ! on output, solution of the linear system
   
contains

!***********************************************************************
!BOP
! !IROUTINE: pcg_verify
! !INTERFACE:

 subroutine pcg_verify(X,B)

! !DESCRIPTION:
! This routine uses a preconditioned conjugate-gradient solver to
! solve the equation $Ax=b$. Both the operator $A$ and preconditioner
! are nine-point stencils. If no preconditioner has been supplied,
! a diagonal preconditioner is applied. Convergence is checked
! every {\em ncheck} steps.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      B ! right hand side of linear system

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(inout) :: &
      X ! on input, an initial guess for the solution
                         ! on output, solution of the linear system

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock ! local block counter

   real (r8) :: &
      rr ! scalar inner product results

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      R   ! residual (b-Ax)

   type (block) :: &
      this_block ! block information for current block

!-----------------------------------------------------------------------
!
! compute initial residual and initialize S
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,this_block)

   do iblock=1,nblocks_tropic
      this_block = get_block(blocks_tropic(iblock),iblock)

      call btrop_operator(R,X,this_block,iblock)
      R(:,:,iblock) = B(:,:,iblock) - R(:,:,iblock)
   end do ! block loop

   !$OMP END PARALLEL DO

   ! ljm tuning
   rr = global_sum(R*R, distrb_tropic, &
                   field_loc_center, RCALCT_B) ! (r,r)

   if (my_task == master_task) then
      if (rr < solv_convrg) then
         write(6,*)'pcg_verify: Success to verify! rr=',rr
      else
         write(6,*)'pcg_verify: Failed to verify! rr=',rr
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine pcg_verify

!***********************************************************************
!BOP
! !IROUTINE: pcg
! !INTERFACE:

 subroutine pcg(X,B)

! !DESCRIPTION:
! This routine uses a preconditioned conjugate-gradient solver to
! solve the equation $Ax=b$. Both the operator $A$ and preconditioner
! are nine-point stencils. If no preconditioner has been supplied,
! a diagonal preconditioner is applied. Convergence is checked
! every {\em ncheck} steps.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      B ! right hand side of linear system

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(inout) :: &
      X ! on input, an initial guess for the solution
                         ! on output, solution of the linear system

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      m, &! local iteration counter
      iblock ! local block counter

   real (r8) :: &
      eta0,eta1,rr ! scalar inner product results

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      R, &! residual (b-Ax)
      S, &! conjugate direction vector
      Q,WORK0,WORK1 ! various cg intermediate results

   character (char_len) :: &
      noconvrg ! error message for no convergence

   type (block) :: &
      this_block ! block information for current block

! ------------------------------------------------ SIMPLE_VERSION ---------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      r_k, p_k, s_k, x_k1, r_k1, p_k1, s_k1, rt_k, st_k, rt_k1, st_k1

   real (r8) :: &
      nu_k, mu_k, a_k, a_k1, a_k2, b_k, b_k1, mu_k1, nu_k1, del_k, gam_k, del_k1, gam_k1


   integer (int_kind) :: k
   
   ! initialize
   r_k = B - simple_A(X)
   rt_k = pcer(r_k)
   nu_k = simple_sum(rt_k * r_k)
   p_k = rt_k
   s_k = simple_A(p_k)
   st_k = pcer(s_k)
   mu_k = simple_sum(p_k * s_k)
   a_k = nu_k / mu_k
   del_k = simple_sum(r_k * st_k)
   gam_k = simple_sum(st_k*s_k)
   a_k1 = 0
   a_k2 = 0
   b_k = 0
   b_k1 = 0
   !

   iter_loop_k: do k = 1, solv_max_iters
      
      !$OMP PARALLEL DO PRIVATE(iblock,this_block)
      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)

         a_k2 = a_k1
         a_k1 = a_k
         b_k1 = b_k
         nu_k1 = nu_k
         del_k1 = del_k
         gam_k1 = gam_k

         x_k1(:,:,iblock) = X(:,:,iblock)
         r_k1(:,:,iblock) = r_k(:,:,iblock)
         rt_k1(:,:,iblock) = rt_k(:,:,iblock)
         p_k1(:,:,iblock) = p_k(:,:,iblock)
         s_k1(:,:,iblock) = s_k(:,:,iblock)
         st_k1(:,:,iblock) = st_k(:,:,iblock)

         X(:,:,iblock) = x_k1(:,:,iblock) + a_k1 * p_k1(:,:,iblock)
         r_k(:,:,iblock) = r_k1(:,:,iblock) - a_k1 * s_k1(:,:,iblock)
         rt_k(:,:,iblock) = rt_k1(:,:,iblock) - a_k1 * st_k1(:,:,iblock)
      end do ! block loop
      !$OMP END PARALLEL DO

      
      ! nu_k = - nu_k1 + (a_k1*a_k1) * gam_k1
      nu_k = nu_k1 - 2 * a_k1 * del_k1 + (a_k1 * a_k1) * gam_k1
      b_k = nu_k / nu_k1

      !$OMP PARALLEL DO PRIVATE(iblock,this_block)
      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)

         p_k(:,:,iblock) = r_k(:,:,iblock) + b_k * p_k1(:,:,iblock)

      end do ! block loop
      !$OMP END PARALLEL DO

      s_k = simple_A(p_k)
      st_k = pcer(s_k)
      

      !$OMP PARALLEL
      !$OMP SECTIONS
      mu_k = simple_sum(p_k * s_k)
      !$OMP END SECTIONS
      !$OMP SECTIONS
      del_k = simple_sum(r_k * st_k)
      !$OMP END SECTIONS
      !$OMP SECTIONS
      gam_k = simple_sum(st_k * s_k)
      !$OMP END SECTIONS
      !$OMP SECTIONS
      nu_k = simple_sum(rt_k * r_k) 
      !$OMP END SECTIONS
      !$OMP END PARALLEL

      a_k = nu_k / mu_k


      ! CHECK
         if (.true.) then

            rr = abs(nu_k)

               ! ljm tuning
            if (my_task == master_task) &
               write(6,*)'  iter k= ',k,' rr= ',rr
            if (rr < solv_convrg) then
               ! ljm tuning
               if (my_task == master_task) &
                  write(6,*)'pcg_iter_loop:iter k= ',k,' rr= ',rr
               solv_sum_iters = k
               exit iter_loop_k
            endif

         endif
      ! ENDCHECK
      

      ! if (my_task == master_task) &
      !    write(6,*) &
      !       'X= ',sum(X), &
      !       'Rk= ',sum(Rk), &
      !       'RHk= ',sum(RHk), &
      !       'Pk= ',sum(Pk), &
      !       'Sk= ',sum(Sk), &
      !       'Vk= ',Vk, &
      !       'Alpha= ',Alphak

   enddo iter_loop_k


   
   ! if (solv_sum_iters == solv_max_iters) then
   !    if (solv_convrg /= c0) then
   !       write(noconvrg,'(a45,i11)') &
   !         'Barotropic solver not converged at time step ', nsteps_total
   !       call exit_POP(sigAbort,noconvrg)
   !    endif
   ! endif

return 
! ------------------------------------------------ SIMPLE_VERSION END---------------------------------------------------------


!-----------------------------------------------------------------------
!
! compute initial residual and initialize S
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,this_block)

   do iblock=1,nblocks_tropic
      this_block = get_block(blocks_tropic(iblock),iblock)

      call btrop_operator(S,X,this_block,iblock)
      R(:,:,iblock) = B(:,:,iblock) - S(:,:,iblock)
      S(:,:,iblock) = c0
   end do ! block loop

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! initialize fields and scalars
!
!-----------------------------------------------------------------------

   call update_ghost_cells(R, bndy_tropic, field_loc_center, &
                                           field_type_scalar)
   eta0 =c1
   solv_sum_iters = solv_max_iters

!-----------------------------------------------------------------------
!
! iterate
!
!-----------------------------------------------------------------------

   iter_loop: do m = 1, solv_max_iters

!-----------------------------------------------------------------------
!
! calculate (PC)r
! diagonal preconditioner if preconditioner not specified
!
!-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock,this_block)

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)

         if (lprecond) then
            call preconditioner(WORK1,R,this_block,iblock)
         else
            where (A0(:,:,iblock) /= c0)
               WORK1(:,:,iblock) = R(:,:,iblock)/A0(:,:,iblock)
            elsewhere
               WORK1(:,:,iblock) = c0
            endwhere
         endif
                                          ! M^{-1} r_i
         WORK0(:,:,iblock) = R(:,:,iblock)*WORK1(:,:,iblock)
         ! r_i^T M^{-1} r_i
      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! update conjugate direction vector s
!
!-----------------------------------------------------------------------

      if (lprecond) &
         call update_ghost_cells(WORK1,bndy_tropic, field_loc_center,&
                                                    field_type_scalar)
      !*** (r,(PC)r)
      ! r_i^T M^{-1} r_i
      eta1 = global_sum(WORK0, distrb_tropic, field_loc_center, RCALCT_B)

      !$OMP PARALLEL DO PRIVATE(iblock,this_block)

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)
                        ! M^{-1} r_i + d (beta_{i+1} / beta_i)
         S(:,:,iblock) = WORK1(:,:,iblock) + S(:,:,iblock)*(eta1/eta0)
        

!-----------------------------------------------------------------------
!
! compute As
!
!-----------------------------------------------------------------------

         call btrop_operator(Q,S,this_block,iblock)
                             ! AS
         WORK0(:,:,iblock) = Q(:,:,iblock)*S(:,:,iblock)
         ! SAS
      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! compute next solution and residual
!
!-----------------------------------------------------------------------

      call update_ghost_cells(Q, bndy_tropic, field_loc_center, &
                                              field_type_scalar)

      eta0 = eta1
      ! r_i^T M^{-1} r_i
      ! r_i^T M^{-1} r_i / SAS
      eta1 = eta0/global_sum(WORK0, distrb_tropic, &
                             field_loc_center, RCALCT_B)
      ! alpha_i
      !$OMP PARALLEL DO PRIVATE(iblock,this_block)

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)
                                             
         X(:,:,iblock) = X(:,:,iblock) + eta1*S(:,:,iblock)
         R(:,:,iblock) = R(:,:,iblock) - eta1*Q(:,:,iblock)

         if (mod(m,solv_ncheck) == 0) then

            call btrop_operator(R,X,this_block,iblock)
            R(:,:,iblock) = B(:,:,iblock) - R(:,:,iblock)
            ! b - Ax
            WORK0(:,:,iblock) = R(:,:,iblock)*R(:,:,iblock)
            ! R^2
         endif
      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! test for convergence
!
!-----------------------------------------------------------------------
      
      
      if (mod(m,solv_ncheck) == 0) then

         call update_ghost_cells(R, bndy_tropic, field_loc_center,&
                                                 field_type_scalar)

         rr = global_sum(WORK0, distrb_tropic, &
                         field_loc_center, RCALCT_B) ! (r,r)

            ! ljm tuning
!            if (my_task == master_task) &
!               write(6,*)'  iter#= ',m,' rr= ',rr
         if (rr < solv_convrg) then
            ! ljm tuning
            if (my_task == master_task) &
               write(6,*)'pcg_iter_loop:iter#=',m,' rr= ',rr
            solv_sum_iters = m
            exit iter_loop
         endif

      endif

   enddo iter_loop

   rms_residual = sqrt(rr*resid_norm)

   
   if (solv_sum_iters == solv_max_iters) then
      if (solv_convrg /= c0) then
         write(noconvrg,'(a45,i11)') &
           'Barotropic solver not converged at time step ', nsteps_total
         call exit_POP(sigAbort,noconvrg)
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine pcg

!***********************************************************************
!BOP
! !IROUTINE: preconditioner
! !INTERFACE:

 subroutine preconditioner(PX,X,this_block,bid)

! !DESCRIPTION:
! This function applies a precomputed preconditioner as a nine-point
! stencil operator.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      X ! array to be operated on

   type (block), intent(in) :: &
      this_block ! block info for this block

   integer (int_kind), intent(in) :: &
      bid ! local block address for this block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(out) :: &
      PX ! nine point operator result

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j ! dummy counters

!-----------------------------------------------------------------------

   PX(:,:,bid) = c0

   do j=this_block%jb,this_block%je
   do i=this_block%ib,this_block%ie
      PX(i,j,bid) = PCNE(i,j,bid)*X(i+1,j+1,bid) + &
                    PCNW(i,j,bid)*X(i-1,j+1,bid) + &
                    PCSE(i,j,bid)*X(i+1,j-1,bid) + &
                    PCSW(i,j,bid)*X(i-1,j-1,bid) + &
                    PCN (i,j,bid)*X(i ,j+1,bid) + &
                    PCS (i,j,bid)*X(i ,j-1,bid) + &
                    PCE (i,j,bid)*X(i+1,j ,bid) + &
                    PCW (i,j,bid)*X(i-1,j ,bid) + &
                    PCC (i,j,bid)*X(i ,j ,bid)
   end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine preconditioner

!***********************************************************************
!BOP
! !IROUTINE: btrop_operator
! !INTERFACE:

function simple_A(X) result (AX)
   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      X ! array to be operated on

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      AX ! nine point operator result (Ax)

   integer (int_kind) :: &
      iblock ! local block counter

   type (block) :: &
      this_block ! block information for current block

   !LINE: 1
   !$OMP PARALLEL DO PRIVATE(iblock,this_block)
      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)
   
        call btrop_operator(AX, X,this_block,iblock)
      end do ! block loop
   !$OMP END PARALLEL DO
   call update_ghost_cells(AX, bndy_tropic, field_loc_center, field_type_scalar)
   !ENDLINE: 1


end function simple_A

function simple_sum(X) result (s)

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      X ! array to be operated on

   real (r8) :: s

   s = global_sum(X, distrb_tropic, &
                         field_loc_center, RCALCT_B) 

end function simple_sum

 subroutine btrop_operator(AX,X,this_block,bid)

! !DESCRIPTION:
! This routine applies the nine-point stencil operator for the
! barotropic solver. It takes advantage of some 9pt weights being
! shifted versions of others.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      X ! array to be operated on

   type (block), intent(in) :: &
      this_block ! block info for this block

   integer (int_kind), intent(in) :: &
      bid ! local block address for this block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(out) :: &
      AX ! nine point operator result (Ax)

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j ! dummy counters

!-----------------------------------------------------------------------

   AX(:,:,bid) = c0

   do j=this_block%jb,this_block%je
   do i=this_block%ib,this_block%ie
      AX(i,j,bid) = A0 (i ,j ,bid)*X(i ,j ,bid) + &
                    AN (i ,j ,bid)*X(i ,j+1,bid) + &
                    AN (i ,j-1,bid)*X(i ,j-1,bid) + &
                    AE (i ,j ,bid)*X(i+1,j ,bid) + &
                    AE (i-1,j ,bid)*X(i-1,j ,bid) + &
                    ANE(i ,j ,bid)*X(i+1,j+1,bid) + &
                    ANE(i ,j-1,bid)*X(i+1,j-1,bid) + &
                    ANE(i-1,j ,bid)*X(i-1,j+1,bid) + &
                    ANE(i-1,j-1,bid)*X(i-1,j-1,bid)
   end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine btrop_operator

include 'matpow.inc.f90'

!***********************************************************************

end module solver_pcg_mod
