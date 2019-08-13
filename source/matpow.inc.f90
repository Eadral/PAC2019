function pcer(X) result (M_X)


    real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
        X, M_X

    M_X = X

end function pcer
