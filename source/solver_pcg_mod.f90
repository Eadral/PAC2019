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
   

   ! interface simple_sum                    ! 虚拟函数名称show
   !    module procedure simple_sum_s       ! 待选择的函数show_int，输入参数为整数
   !    module procedure simple_sum_v ! 待选择的函数show_character 输入参数为字符
   ! end interface

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

   ! --------------------------------------------- MINE ---------------------------------------------------------
   
   integer (int_kind), parameter :: step = 1
   
   real (r8), dimension(nx_block,ny_block,max_blocks_tropic, 0 : step+1) :: &
      Vk

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic, step) :: &
      Rk1s, Rk1s_1
   
   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      Xi, Ri, Ri1, Vsk1, Rsk1, &
      T, &
      Xj1, Xj_1, Wj, Rj, Rj1, Rj_1

   real (r8), dimension(ny_block,max_blocks_tropic) :: &
      L, Uj, Vj, Gammaj_1, Gammaj, Rouj, Rouj_1, Rouk, Gammak, Uj_1

   
   integer (int_kind) :: k, j, rs, vs

   real (r8) :: &
      Ts1, Ts2


   real (r8), dimension(ny_block,max_blocks_tropic, 2*step+1, 1) :: &
      Dj, Dj_1, Dj_2
   
   real (r8), dimension(step, step) :: &
      Tk_1, Tk


   real (r8), dimension(ny_block,max_blocks_tropic, step) :: Rouk1s, Gammak1s

   real (r8), dimension(step+1, step) :: &
      Bk
   ! init

   !LINE: 1
   Xj_1 = c0
   Rj_1 = c0

   Xj1 = X
   Rj1 = B - simple_A(X)

   ! for k = 0
   Rouk1s = 1
   Gammak1s = 1

   Rouk = 1
   Gammak = 1
   
   Tk_1 = 0
   Rk1s_1 = c0



   capcg_loop_k: do k = 0, solv_max_iters

      ! may no need
      ! Rj1 = B - simple_A(X)

      Vsk1 = Rj1

      Vk = matpow(Vsk1, step)
      Bk = compute_B(step)

      ! for j = 2: Dj_1, Dj_2
      Dj = compute_d(step, Tk_1, Bk, k, Rouk, Gammak, 1, Rouj_1, Gammaj_1, Dj_1, Dj_2)
      Dj_1 = compute_d(step, Tk_1, Bk, k, Rouk, Gammak, 0, Rouj_1, Gammaj_1, Dj_1, Dj_2)
 
      capcg_loop_j: do j = 1, step
         X = Xj1
         Rj = Rj1

         ! if (my_task == master_task) &
         !    write(6,*)'  iter k= ',k,'j=',j,'Dj_1= ',Dj_1
         ! if (my_task == master_task) &
         !    write(6,*)'  iter k= ',k,'j=',j,'Dj_2= ',Dj_2
         Dj = compute_d(step, Tk_1, Bk, k, Rouk, Gammak, j, Rouj_1, Gammaj_1, Dj_1, Dj_2)

         ! if (my_task == master_task) &
         !    write(6,*) Dj

         Wj = c0
         do rs = 1, step
            Wj = Wj + Rk1s_1(:,:,:, rs) * broad_rv(Dj(:,:, rs, 1))
         enddo 

         ! if (my_task == master_task) &
         !    write(6,*)'  iter k= ',k,'j=',j,'Wj= ',sum(Wj)

         do vs = 1, step+1
            Wj = Wj + Vk(:,:,:, vs) * broad_rv(Dj(:,:, step + vs, 1))

            ! if (my_task == master_task) &
            ! write(6,*)'  iter k= ',k,'j=',j,'Dj= ',Dj(step + vs, 1)
         enddo
         
         ! Wj = simple_A(Rj)
         
         ! rr = simple_sum(Wj-simple_A(Rj))
         ! if (my_task == master_task) &
         ! write(6,*)'diff=', rr

         !LINE: 5
         Uj = simple_sum_v(Rj*Rj)

         ! if (my_task == master_task) &
         ! write(6,*)'uj=', Uj

         ! CHECK
         if (.true.) then

            rr = simple_sum_s(Rj*Rj)

               ! ljm tuning
            if (my_task == master_task) &
               write(6,*)'  iter k= ',k,'j=',j,' rr= ',rr
            if (rr < solv_convrg) then
               ! ljm tuning
               if (my_task == master_task) &
                  write(6,*)'pcg_iter_loop:iter k= ',k,' rr= ',rr
               solv_sum_iters = k
               exit capcg_loop_k
            endif

         endif
         ! ENDCHECK
         
         !LINE: 6
         Vj = simple_sum_v(Wj*Rj)

         !LINE: 7
         Gammaj = Uj / Vj

         !LINE: 8
         if ( k == 0 .and. j == 1 ) then
            Rouj = c1
         else
            Rouj = c1 / (c1 - ( (Gammaj/Gammaj_1) * (Uj/Uj_1) * (1/Rouj_1) ) )
         end if
      
         !LINE: 12
         Xj1 = broad_rv(Rouj) * (X + broad_rv(Gammaj)*Rj) + (1 - broad_rv(Rouj))*Xj_1
         !LINE: 13
         Rj1 = broad_rv(Rouj) * (Rj - broad_rv(Gammaj)*Wj) + (1 - broad_rv(Rouj))*Rj_1
         
         ! update sk+j
         Rouk1s(:,:, j) = Rouj
         Gammak1s(:,:, j) = Gammaj
         Rk1s(:,:,:, j) = Rj

         ! update j_1
         Gammaj_1 = Gammaj
         Uj_1 = Uj
         Rouj_1 = Rouj
         Xj_1 = X
         Rj_1 = Rj

         Dj_2 = Dj_1
         Dj_1 = Dj
         
      enddo capcg_loop_j

      ! update k
      Rouk = Rouj
      Gammak = Gammaj


      Tk =  compute_T(step, Rouk1s, Gammak1s)
      ! update k_1
      Tk_1 = Tk
      Rk1s_1 = Rk1s
   enddo capcg_loop_k

   ! if (solv_sum_iters == solv_max_iters) then
   !    if (solv_convrg /= c0) then
   !       write(noconvrg,'(a45,i11)') &
   !         'Barotropic solver not converged at time step ', nsteps_total
   !       call exit_POP(sigAbort,noconvrg)
   !    endif
   ! endif
   return
! ------------------------------------------------ MINE END---------------------------------------------------------


! ------------------------------------------------ SIMPLE_VERSION ---------------------------------------------------------


   R = B - simple_A(X)
   S = c0

   eta0 =c1
   solv_sum_iters = solv_max_iters

   iter_loop_m: do m = 1, solv_max_iters


      where (A0 /= c0)
         WORK1 = R/A0
      elsewhere
         WORK1 = c0
      endwhere
 
                                 ! M^{-1} r_i
      WORK0 = R*WORK1
         ! r_i^T M^{-1} r_i

      !*** (r,(PC)r)
      ! r_i^T M^{-1} r_i
      eta1 = simple_sum_s(WORK0)

      S = WORK1 + S*(eta1/eta0)
      Q = simple_A(S)
                             ! AS
      WORK0 = Q*S
      ! SAS

      eta0 = eta1
      ! r_i^T M^{-1} r_i
      ! r_i^T M^{-1} r_i / SAS
      eta1 = eta0 / simple_sum_s(WORK0)
      ! alpha_i
    
                                             
      X = X + eta1*S
      R = R - eta1*Q

      if (mod(m,solv_ncheck) == 0) then

         R = simple_A(X)
         R = B - R
         ! b - Ax
         WORK0 = R*R
         ! R^2

         rr = simple_sum_s(WORK0)

            ! ljm tuning
!            if (my_task == master_task) &
!               write(6,*)'  iter#= ',m,' rr= ',rr
         if (rr < solv_convrg) then
            ! ljm tuning
            if (my_task == master_task) &
               write(6,*)'pcg_iter_loop:iter#=',m,' rr= ',rr
            solv_sum_iters = m
            exit iter_loop_m
         endif

      endif

   enddo iter_loop_m

   rms_residual = sqrt(rr*resid_norm)

   
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

function simple_sum_s(X) result (s)

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      X ! array to be operated on

   real (r8) :: s

   s = global_sum(X, distrb_tropic, &
                         field_loc_center, RCALCT_B) 

end function simple_sum_s

function simple_sum_v(X) result (V)

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      X ! array to be operated on

   real (r8), dimension(ny_block, max_blocks_tropic) &
      :: V

   integer (int_kind) :: inx
   
   V = 0
   do inx = 1, nx_block
      V = V + X(inx, :, :)
   enddo

end function simple_sum_v

! function simple_sum_rv_s(X) result (s)

!    real (r8), dimension(ny_block,max_blocks_tropic), &
!       intent(in) :: &
!       X ! array to be operated on

!       real (r8) :: s

!       s = global_sum(X, distrb_tropic, &
!                             field_loc_center, RCALCT_B) 
   
! end function simple_sum_rv_s

function broad_rv(X) result (V)

   real (r8), dimension(ny_block,max_blocks_tropic), &
      intent(in) :: &
      X ! array to be operated on

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: V

     V = spread(X, 1, nx_block)
   
end function broad_rv

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
