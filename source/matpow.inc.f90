subroutine matpow(Vk, V1, step)

    integer (int_kind) :: step

    real (r8), dimension(nx_block,ny_block,max_blocks_tropic), intent(in) :: &
        V1

    real (r8), dimension(nx_block,ny_block,max_blocks_tropic, 0 : step+1), intent(out) :: &
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

end subroutine matpow