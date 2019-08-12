function matpow(V1, step) result (Vk)

    integer (int_kind) :: step

    real (r8), dimension(nx_block,ny_block,max_blocks_tropic), intent(in) :: &
        V1

    real (r8), dimension(nx_block,ny_block,max_blocks_tropic, 0 : step+1) :: &
        Vk
    
    integer (int_kind) :: k
    ! real (r8), dimension(nx_block,ny_block,max_blocks_tropic):: &
    !     T

    Vk(:,:,:, 0) = c0
    Vk(:,:,:, 1) = V1

    do k = 1, step
        Vk(:,:,:, k+1) = simple_A(Vk(:,:,:, k))
    enddo

end function matpow

function compute_d(s, Tk_1, Bk, k, Rouk, Gammak, j, Rouj_1, Gammaj_1, Dj_1, Dj_2) result (D)

    integer (int_kind) :: s, k, j, iblock, iny

    real (r8), dimension(ny_block,max_blocks_tropic, 2*s+1, 1) :: &
        D, DT1, DT2, Dj_1, Dj_2

    real (r8), dimension(s, s) :: &
        Tk_1
        
    real (r8), dimension(s+1, s) :: &
        Bk

    real (r8), dimension(ny_block,max_blocks_tropic) :: Rouk, GammaK, Rouj_1, Gammaj_1



    if ( k == 0 .and. j == 0 ) then
        D = c0
    else if ( k > 0 .and. j == 0) then
        D = c0
        D(:,:, s, 1) = c1
    else if (j == 1) then
        D = c0
        D(:,:,s+1, 1) = Bk(1, 1)
        D(:,:,s+2, 1) = Bk(2, 1)
    else

        ! if (my_task == master_task) &
        ! write(6,*)'  iter k= ',k,'j=',j,'Tk_1= ',sum(Tk_1)
        ! if (my_task == master_task) &
        ! write(6,*)'  iter k= ',k,'j=',j,'Dj_1= ',sum(Dj_1)
     

        DT1 = 0
        do iblock=1,nblocks_tropic
            do iny=1, ny_block
                DT1(iny,iblock, 1:s, 1) = matmul(Tk_1,  Dj_1(iny,iblock, 1:s, 1) )
            enddo
        enddo
        
        DT1(:,:, s+1, 1) = - (Dj_1(:,:,s, 1) / (Rouk * GammaK))
        
        ! if (my_task == master_task) &
        ! write(6,*)'  iter k= ',k,'j=',j,'dt1= ',sum(DT1)

        DT2 = 0
        do iblock=1,nblocks_tropic
            do iny=1, ny_block
                DT2(iny,iblock, s+1:2*s+1, 1) = matmul(Bk, Dj_1(iny,iblock,s+1:2*s, 1) )
            enddo
        enddo

        ! if (my_task == master_task) &
        ! write(6,*)'  iter k= ',k,'j=',j,'dt2= ',sum(DT2)
        
        do iblock=1,nblocks_tropic
            do iny=1, ny_block
                D(iny,iblock,:,1) = Rouj_1(iny,iblock)*Dj_1(iny,iblock,:,1) + (1 - Rouj_1(iny,iblock))*Dj_2(iny,iblock,:,1) - &
                    Rouj_1(iny,iblock)*Gammaj_1(iny,iblock)*DT1(iny,iblock,:,1) + &
                    Rouj_1(iny,iblock)*Gammaj_1(iny,iblock)*DT2(iny,iblock,:,1)
            enddo
        enddo
        

        ! if (my_task == master_task) &
        ! write(6,*)'  iter k= ',k,'j=',j,'D= ',sum(D)
        

    end if

end function compute_d

function compute_T(s, Rouk1s, Gammak1s) result (T)

    integer (int_kind) :: i, s

    real (r8), dimension(s, s) :: T

    real (r8), dimension(s+1, s) :: Tu

    real (r8), dimension(s), intent(in) :: Rouk1s, Gammak1s 

    real (r8), dimension(s+1, s) :: Tt

    real (r8), dimension(s, s) :: diag


    Tt = 0
    diag = 0

    do i = 1, s
        diag(i, i) = 1 / (Rouk1s(i)*Gammak1s(i))
    enddo

    Tt(1, 1) = Rouk1s(1)
    Tt(2, 1) = -1
    do i = 2, s
        Tt(i-1, i) = 1 - Rouk1s(i)
        Tt(i, i) = Rouk1s(i)
        Tt(i+1, i) = -1
    enddo

    Tu = matmul(Tt, diag)

    T = Tu(1:s, 1:s)

end function compute_T

function compute_B(s) result(B)

    integer (int_kind) :: s

    real (r8), dimension(s+1, s) :: &
        B

    ! real (r8), dimension(nx_block,ny_block,max_blocks_tropic), intent(in) :: &
       
    integer (int_kind) :: k

    B = c0
    do k = 1, s
        B(k+1, k) = 1
        ! B(k, k) = 0.000000001
    enddo
end function compute_B