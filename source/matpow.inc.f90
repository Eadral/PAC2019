function matpow(V1, step) result (Vk)

    integer (int_kind) :: step

    real (r8), dimension(nx_block,ny_block,max_blocks_tropic), intent(in) :: &
        V1

    real (r8), dimension(nx_block,ny_block,max_blocks_tropic, 0 : step+1) :: &
        Vk
    
    integer (int_kind) :: k
    real (r8), dimension(nx_block,ny_block,max_blocks_tropic):: &
        T

    Vk(:,:,:, 0) = c0
    Vk(:,:,:, 1) = V1

    do k = 1, step
        T = simple_A(Vk(:,:,:, k))
        Vk(:,:,:, k+1) = T 
    enddo

end function matpow

! function compute_d(Bk, k, Rouk, Gammak, j, Rouj_1, Gammaj_1, Dj_1, Dj_2) result (D)

!     real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
!         D

!     real (r8), dimension(nx_block,ny_block,max_blocks_tropic), intent(in) :: &
        

!     real (r8), intent(in) :: 

!     if ( k == 0 .and. j == 0 ) then
!         D = c0
!     else if ( k > 0 .and j == 0) then
!         stop
!     else if (j == 1) then
!         D = c0

!     else

!     end if

! end function compute_d


! function compute_B(step) result(B)

!     real (r8), dimension(step+1, step) :: &
!         B

!     ! real (r8), dimension(nx_block,ny_block,max_blocks_tropic), intent(in) :: &
       
!     integer (int_kind) :: k

!     B = c0
!     do k = 1, step
!         B(k+1, k) = 1
!     enddo
! end function compute_B