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

   integer (int_kind), parameter :: step = 3

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
     r0, x0, p0, AX

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
     ri, xi, pi


   integer (int_kind) :: &
      k, its, gi, gj, j, inner, outer, jj, i

   !
   real (r8), dimension(step, 1) :: &
      alp, gam

   real (r8), dimension(step-1, 1) :: &
      bet

   real (r8), dimension(step+1, step) :: &
      T

   
   !

   !
   real (r8), dimension(nx_block,ny_block,max_blocks_tropic,step+1) :: &
      P
   
   real (r8), dimension(nx_block,ny_block,max_blocks_tropic,step) :: &
      Rs
   !

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic,2*step+1) :: &
      PR
   

   real (r8), dimension(2*step+1, 2*step+1) :: &
      Tt

   real (r8), dimension(2*step+1, step+1) :: &
      p_c, r_c, x_c

   real (r8), dimension(2*step+1, 2*step+1) :: &
      G

   real (r8), dimension(2*step*step+2*step+2) :: &
      Gl

   real (r8), dimension(solv_max_iters) :: &
      alpha, beta
   
   real (r8), dimension(step) :: &
      D

   real (r8) :: &
      temp1, ts, te

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      P_P, P_R, R_R

   integer (int_kind) :: kc
   
   !$OMP PARALLEL DO PRIVATE(iblock,this_block)  
   do iblock=1,nblocks_tropic
      this_block = get_block(blocks_tropic(iblock),iblock)
      
      call btrop_operator(AX, X, this_block,iblock)

   enddo
   !$OMP END PARALLEL DO

   call update_ghost_cells(AX, bndy_tropic, field_loc_center, field_type_scalar)

   !$OMP PARALLEL DO PRIVATE(iblock,this_block)  
   do iblock=1,nblocks_tropic
      this_block = get_block(blocks_tropic(iblock),iblock)
      

      ri(:,:,iblock) = B(:,:,iblock) - AX(:,:,iblock)
      pi(:,:,iblock) = ri(:,:,iblock)
   enddo
   !$OMP END PARALLEL DO
   ! AX = simple_A(X)

   ! r0 = B - AX
   ! p0 = r0
   ! xi(:,:,:) = X
   ! ri(:,:,:) = B - AX
   ! pi(:,:,:) = ri(:,:,:)

   k = 0

   its = 0

   call basisparams(step, alp, bet, gam, T)

   ! solv_max_iters = 180

   iters: do while (its < solv_max_iters)

      ! call cpu_time(ts)


      ! PR(:,:,:, 1:step+1) = computeBasis(pi(:,:,:), step, alp, bet, gam)

      ! PR(:,:,:, step+2:2*step+1) = computeBasis(ri(:,:,:), step-1, alp, bet, gam)

         !$OMP PARALLEL DO PRIVATE(iblock,this_block)  
         do iblock=1,nblocks_tropic
            this_block = get_block(blocks_tropic(iblock),iblock)
            
            PR(:,:,iblock, 1) = pi(:,:,iblock)

            PR(:,:,iblock, step+2) = ri(:,:,iblock)
         enddo
         !$OMP END PARALLEL DO

      

         do jj = 1, step-1
            !$OMP PARALLEL DO PRIVATE(iblock,this_block)    
            do iblock=1,nblocks_tropic
               this_block = get_block(blocks_tropic(iblock),iblock)
      
               PR(:,:,iblock, jj+1) = c0
               PR(:,:,iblock, step+1+jj+1) = c0
               
               !$OMP PARALLEL DO
               do j=this_block%jb,this_block%je
               !DIR$ IVDEP UNROLL
                  !$OMP PARALLEL DO
               do i=this_block%ib,this_block%ie
                  PR(i,j,iblock, jj+1) = A0 (i ,j ,iblock)*PR(i ,j ,iblock,jj) + &
                                AN (i ,j ,iblock)*PR(i ,j+1,iblock,jj) + &
                                AN (i ,j-1,iblock)*PR(i ,j-1,iblock,jj) + &
                                AE (i ,j ,iblock)*PR(i+1,j ,iblock,jj) + &
                                AE (i-1,j ,iblock)*PR(i-1,j ,iblock,jj) + &
                                ANE(i ,j ,iblock)*PR(i+1,j+1,iblock,jj) + &
                                ANE(i ,j-1,iblock)*PR(i+1,j-1,iblock,jj) + &
                                ANE(i-1,j ,iblock)*PR(i-1,j+1,iblock,jj) + &
                                ANE(i-1,j-1,iblock)*PR(i-1,j-1,iblock,jj)

                  PR(i,j,iblock, step+1+jj+1) = A0 (i ,j ,iblock)*PR(i ,j ,iblock, step+1+jj) + &
                                AN (i ,j ,iblock)*PR(i ,j+1,iblock, step+1+jj) + &
                                AN (i ,j-1,iblock)*PR(i ,j-1,iblock, step+1+jj) + &
                                AE (i ,j ,iblock)*PR(i+1,j ,iblock, step+1+jj) + &
                                AE (i-1,j ,iblock)*PR(i-1,j ,iblock, step+1+jj) + &
                                ANE(i ,j ,iblock)*PR(i+1,j+1,iblock, step+1+jj) + &
                                ANE(i ,j-1,iblock)*PR(i+1,j-1,iblock, step+1+jj) + &
                                ANE(i-1,j ,iblock)*PR(i-1,j+1,iblock, step+1+jj) + &
                                ANE(i-1,j-1,iblock)*PR(i-1,j-1,iblock, step+1+jj)
               end do
               !$OMP END PARALLEL DO
               end do
               !$OMP END PARALLEL DO
            
            ! call btrop_operator(PR(:,:,:, jj+1), PR(:,:,:, jj), this_block,iblock)
            ! call btrop_operator(PR(:,:,:, step+1+jj+1), PR(:,:,:, step+1+jj), this_block,iblock)
            enddo
            !$OMP END PARALLEL DO
            call update_ghost_cells(PR(:,:,:, jj+1), bndy_tropic, field_loc_center, field_type_scalar)
            call update_ghost_cells(PR(:,:,:, step+1+jj+1), bndy_tropic, field_loc_center, field_type_scalar)
         enddo

         jj = step
         !$OMP PARALLEL DO PRIVATE(iblock,this_block)  
         do iblock=1,nblocks_tropic
            this_block = get_block(blocks_tropic(iblock),iblock)
   
         call btrop_operator(PR(:,:,:, jj+1), PR(:,:,:, jj), this_block,iblock)
         enddo
         !$OMP END PARALLEL DO
         call update_ghost_cells(PR(:,:,:, jj+1), bndy_tropic, field_loc_center, field_type_scalar)

      
      Tt = 0
      Tt(1:step+1, 1:step) = T
      Tt(step+2:2*step+1, step+2:2*step) = T(1:step, 1:(step-1))

      p_c = 0
      p_c(1, 1) = 1

      r_c = 0
      r_c(step+2, 1) = 1
      
      x_c = 0
      
      ! PR(:,:,:, 1:step+1) = P
      ! PR(:,:,:, step+2:2*step+1) = Rs
      ! G = 0
      ! do gi = 1, 2*step+1
      !    do gj = 1, 2*step+1
      !       G(gi, gj) = simple_sum( PR(:,:,:, gi) * PR(:,:,:, gj) )
      !    enddo
      ! enddo
      ! call cpu_time(te)
      ! if (my_task == master_task) write(6, *)'3', (te-ts)*1000
      ! call cpu_time(ts)



      G = 0

      !$OMP PARALLEL DO PRIVATE(iblock,this_block) 

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)

         do gi = 1, 2*step+1
            do gj = gi, 2*step+1
              
               P_R(:,:,iblock) = PR(:,:,iblock, gi) * PR(:,:,iblock, gj)
               temp1 = local_sum( P_R, iblock )
               G(gi, gj) = G(gi, gj) + temp1
               if (.not.(gi == gj)) &
               G(gj, gi) = G(gj, gi) + temp1
            enddo
         enddo

      end do ! block loop

      ! call cpu_time(te)
      ! if (my_task == master_task) write(6, *)'4', (te-ts)*1000
      ! call cpu_time(ts)
      G = global_sum(G, distrb_tropic)
      ! call cpu_time(te)
      ! if (my_task == master_task) write(6, *)'5', (te-ts)*1000
      ! call cpu_time(ts)

      ! stop

      do j = 1, step
         
         if (its >= solv_max_iters) then
            exit iters
         endif 

         its = its + 1

         D(j) = gram(r_c(:,j), G, r_c(:,j), step) 

         alpha(its) = D(j) / gram_Tt(p_c(:,j), G, Tt, p_c(:,j), step )  

         x_c(:,j+1) = x_c(:,j) + alpha(its)*p_c(:,j)

         r_c(:,j+1) = r_c(:,j) - matmul(alpha(its)*Tt, p_c(:,j) )

         beta(its) =  gram(r_c(:,j+1), G, r_c(:,j+1), step) /  D(j)

         p_c(:,j+1) = r_c(:,j+1) + beta(its)*p_c(:,j)


      enddo

      ! call cpu_time(te)
      ! if (my_task == master_task) write(6, *)'6', (te-ts)*1000
      ! call cpu_time(ts)

      !$OMP PARALLEL DO PRIVATE(iblock,this_block,P_R)

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)

         ri(:,:,iblock) = 0

         pi(:,:,iblock) = 0

         do kc = 1, 2*step+1

            X(:,:,iblock) = X(:,:,iblock) + x_c(kc,step+1) * PR(:,:,iblock,kc)

            ri(:,:,iblock) = ri(:,:,iblock) + r_c(kc,step+1) * PR(:,:,iblock,kc)

            pi(:,:,iblock) = pi(:,:,iblock) +  p_c(kc,step+1) * PR(:,:,iblock,kc)

         enddo

         ! X(:,:,iblock) = xi(:,:,iblock)

      end do ! block loop

      !$OMP END PARALLEL DO


      k = k+1

      
       ! CHECK
      if (mod(its, solv_ncheck) == 0) then

         rr = simple_sum( ri )
         rr = rr * rr / 10

            ! ljm tuning
         if (my_task == master_task) &
            write(6,*)'  iter its= ',its,' rr= ',rr
         if (rr < solv_convrg) then
            ! ljm tuning
            if (my_task == master_task) &
               write(6,*)'pcg_iter_loop:iter its= ',its,' rr= ',rr
            solv_sum_iters = its
            exit iters
         endif

      endif
   ! ENDCHECK

      ! call cpu_time(te)
      ! if (my_task == master_task) write(6, *)'7', (te-ts)*1000

   enddo iters


return 
end subroutine pcg
! ------------------------------------------------ SIMPLE_VERSION END---------------------------------------------------------

subroutine basisparams(s, alp, bet, gam, T)

   integer (int_kind) :: &
      k, s

   real (r8), dimension(s, 1), intent(out):: &
      alp, gam

   real (r8), dimension(s-1, 1), intent(out) :: &
      bet

   real (r8), dimension(s+1, s), intent(out) :: &
      T
   
   ! monomial
   alp = 0 
   bet = 0
   gam = 1

   T = 0
   do k = 1, s
      T(k+1, k) = gam(k, 1)
   enddo


end subroutine basisparams

function computeBasis(x, s, alp, bet, gam) result (V)

   integer (int_kind) :: &
      j, s, iblock

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic,s+1) :: &
     V

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
     x

   real (r8), dimension(s, 1) :: &
      alp, gam

   real (r8), dimension(s-1, 1) :: &
      bet

   real (r8) :: ts, te, t1, t2

   type (block) :: &
      this_block ! block information for current block

   V(:,:,:, 1) = x

!   t1 = 0
!   t2 = 0

   if (s > 0) then

   
      ! V(:,:,:, 2) = (1/gam(1, 1)) * simple_A( V(:,:,:, 1) )
      
      do j = 1, s

      !        !LINE: 1
      !$OMP PARALLEL DO PRIVATE(iblock,this_block)  
      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)

         ! call cpu_time(ts)

      call btrop_operator(V(:,:,:, j+1), V(:,:,:, j), this_block,iblock)
      ! call cpu_time(te)

      ! t1 = t1 + te-  ts
      end do ! block loop
      !$OMP END PARALLEL DO
      ! call cpu_time(ts)
      call update_ghost_cells(V(:,:,:, j+1), bndy_tropic, field_loc_center, field_type_scalar)
      ! call cpu_time(te)

      ! t2 = t2 + te - ts
      !ENDLINE: 1


         ! V(:,:,:, j+1) =  simple_A( V(:,:,:, j) )
      enddo

      ! if (my_task == master_task) write(6, *)'3.1', (t1)*1000

      ! if (my_task == master_task) write(6, *)'3.2', (t2)*1000
   endif

end function computeBasis

function eye() result (Y)

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      Y

   where (A0 /= c0)
      Y = 1
   elsewhere
      Y = c0
   endwhere

endfunction eye

function gram(a, G, b, s) result (r)

   integer (int_kind) :: &
      s, ki, kj

   real (r8), dimension(2*s+1) :: &
      a, b

   real (r8), dimension(2*s+1, 2*s+1) :: &
      G
   
   real (r8) :: &
      r, r1, r2

   real (r8), dimension(1, 2*s+1) :: &
      aa

   real (r8), dimension(2*s+1, 1) :: &
      bb

   real (r8), dimension(1, 1) :: &
      rr

   aa = spread(a, dim=1, ncopies=1)

   bb = spread(b, dim=2, ncopies=1)

   rr = matmul( matmul( aa, G ), bb )

   r = rr(1, 1)


   ! r2 = 0
   ! do ki = 1, 2*s+1
   !    r1 = 0
   !    do kj = 1, 2*s+1
   !       r1 = r1 + a(kj) * G(kj, ki)
   !    enddo
   !    r2 = r2 + r1 * b(ki)
   ! enddo

   ! r = r2

endfunction gram

function gram_Tt(a, G, Tt, b, s) result (r)

   integer (int_kind) :: &
      s, ki, kj, kk

   real (r8), dimension(2*s+1) :: &
      a, b

   real (r8), dimension(2*s+1, 2*s+1) :: &
      G, Tt
   
   real (r8) :: &
      r

   real (r8), dimension(1, 2*s+1) :: &
      aa

   real (r8), dimension(2*s+1, 1) :: &
      bb

   real (r8), dimension(1, 1) :: &
      rr

   aa = spread(a, dim=1, ncopies=1)

   bb = spread(b, dim=2, ncopies=1)

   ! rr = matmul( matmul( aa, G ), bb )
   rr = matmul(matmul(matmul(aa, G), Tt), bb )

   r = rr(1, 1)

endfunction gram_Tt

function row_v(x) result (y)

   real (r8), dimension(:) :: x

   real (r8), dimension(1, size(x, dim=1)) :: y

   y = spread(x, dim=1, ncopies=1)

endfunction row_v

function col_v(x) result (y)

   real (r8), dimension(:) :: x

   real (r8), dimension(size(x, dim=1), 1) :: y

   y = spread(x, dim=2, ncopies=1)

end function col_v

function coeff_add(c, X) result (V)

   real (r8), dimension(:) :: c

   real (r8), dimension(nx_block,ny_block, size(c, dim=1)) :: &
      X

   real (r8), dimension(nx_block,ny_block) :: &
      V

   integer (int_kind) :: k

   V = 0
   do k = 1,  size(c, dim=1)
      V = V + c(k) * X(:,:, k)
   enddo

   

end function coeff_add

!***********************************************************************
!BOP
! !IROUTINE: btrop_operator
! !INTERFACE:

! ORIGIN

function local_sum(X, iblock)

   real (r8), dimension(:,:,:), intent(in) :: &
      X                    ! array to be summed

   type (distrb) :: &
      dist                 ! block distribution for array X

   integer (int_kind) :: &
      i,j,n,             &! local counters
      ib,ie,jb,je,       &! beg,end of physical domain
      iblock               ! block location

   real (r8) ::          &
      local_sum           ! sum of all local blocks

   dist = distrb_tropic

   
   ! do n=1,dist%local_block_num
      ! iblock = n
      call get_block_parameter(dist%local_block_ids(iblock),ib=ib,ie=ie,jb=jb,je=je)
         do j=jb,je
         do i=ib,ie
            local_sum = local_sum + X(i,j,iblock)
         end do
         end do
   ! end do !block loop

endfunction local_sum

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

 subroutine btrop_operator(AX,X,this_block,iblock)

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
      iblock ! local block address for this block

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

   AX(:,:,iblock) = c0

   do j=this_block%jb,this_block%je
   do i=this_block%ib,this_block%ie
      AX(i,j,iblock) = A0 (i ,j ,iblock)*X(i ,j ,iblock) + &
                    AN (i ,j ,iblock)*X(i ,j+1,iblock) + &
                    AN (i ,j-1,iblock)*X(i ,j-1,iblock) + &
                    AE (i ,j ,iblock)*X(i+1,j ,iblock) + &
                    AE (i-1,j ,iblock)*X(i-1,j ,iblock) + &
                    ANE(i ,j ,iblock)*X(i+1,j+1,iblock) + &
                    ANE(i ,j-1,iblock)*X(i+1,j-1,iblock) + &
                    ANE(i-1,j ,iblock)*X(i-1,j+1,iblock) + &
                    ANE(i-1,j-1,iblock)*X(i-1,j-1,iblock)
   end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine btrop_operator

include 'matpow.inc.f90'

!***********************************************************************

end module solver_pcg_mod
