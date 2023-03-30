! Program to assembly the stiffness matrix associated 
! with the elliptic PDE
!
!   a u_xx + 2b u_xy +c u_yy = A +Bx^2 + Cy^2 + Exy + D(xy)^2   in [-l,1]x[-h,h]
!   u = G in fr([-l,1]x[-h,h])

program smatrix
    implicit none
    integer :: i, j, k, n_x, n_y, p
    real(8) :: a, b, c, h_x, h_y, A_x, B_y, C_xy, D_xy, l, h
    real(8), allocatable :: S(:,:)
    
    n_x = 4 
    n_y = 4
    l = 0.5d0
    h = 0.5d0
    h_x = 2d0*l/real(n_x, kind=8)
    h_y = 2d0*h/real(n_y, kind=8)
    
    A_x = a/(h_x)**2
    B_y = 2d0*b/(h_y)**2
    C_xy = c/(4d0*h_x*h_y)
    D_xy = -A_x - B_y
    
    ! Generating the stiffness matrix S :

    p = (n_x-1)*(n_y-1)
    allocate(S(p, p))
    S = 0d0
    

    do j = 1, n_y-1
        do i = 1, n_x-1
            k = (n_x-1)*(j-1) + i
            if ( i > 1 ) then
                S(k, k-1) = A_x
                S(k, n_x + k-2) = -C_xy
            end if
            S(k,k) = D_xy
            S(k, n_x + k) = C_xy
            S(k,k+1) = A_x
            S(k, n_x + k-1) = B_y
        end do
    end do
    print *, S(2,:)
    deallocate(S)
end program smatrix