program project2
    implicit none
    integer :: n_x, n_y
    real(8) :: a, b, c, l, h, T, Q, M, E, D, G

    interface
        subroutine generate_matrices(n_x,n_y,h,l,a,b,c,T,Q,M,E,D,G)
            integer :: i, j, k, p
            integer, intent(in) :: n_x, n_y
            real(8), intent(in) :: a, b, c, h, l, T, Q, M, E, D, G
            real(8) :: h_x, h_y, A_x, B_y, C_xy, D_xy
            real(8), allocatable :: S(:,:), x(:), y(:), F(:,:)
        end subroutine generate_matrices
    end interface

    ! Getting data from user

    print *, 'n_x e n_y '
    read *, n_x, n_y
    print *, 'h e l  '
    read *, h, l
    print *, 'a, b e c  '
    read *, a, b, c
    print *, 'T, Q, M, E, D e G '
    read *, T, Q, M, E, D, G

    call generate_matrices(n_x,n_y,h,l,a,b,c,T,Q,M,E,D,G)

end program project2

subroutine generate_matrices(n_x,n_y,h,l,a,b,c,T,Q,M,E,D,G)
    integer :: i, j, k, p
    integer, intent(in) :: n_x, n_y
    real(8), intent(in) :: a, b, c, h, l, T, Q, M, E, D, G
    real(8) :: h_x, h_y, A_x, B_y, C_xy, D_xy
    real(8), allocatable :: S(:,:), x(:), y(:), F(:,:)
    
    h_x = 2d0*l/real(n_x, kind=8)
    h_y = 2d0*h/real(n_y, kind=8)
    p = (n_x-1)*(n_y-1)
    A_x = a/(h_x)**2
    B_y = 2d0*b/(h_y)**2
    C_xy = c/(4d0*h_x*h_y)
    D_xy = -A_x - B_y

    allocate(x(n_x-1))
    allocate(y(n_y-1))

    ! Generating interior mesh poitns:

    do i = 1, n_x-1
        x(i) = -h + h_x*real(i, kind=8) 
    end do

    do j = 1, n_y-1
        y(j) = -l + h_y*real(j, kind=8) 
    end do

    call generate_actions_matrix()
    call generate_stiffness_matrix()
    call matrices_to_csv()
    
    contains
    subroutine generate_actions_matrix()
        allocate(F(p,1))
    
        do j = 1, n_y-1
            do i = 1, n_x-1
                k = (n_x-1)*(j-1) + i
                F(k,1) = T + Q*x(i)**2 + M*y(j)**2 + E*x(i)*y(j) + D*(x(i)*y(i))**2
            end do
        end do
        
        F(1,1) = F(1,1) + (-A_x-B_y+C_xy)*G     ! first entry of F is different 
    
        do i = 1, n_x-1             ! here j = 1
            k = i
            F(k,1) = F(k,1) - B_y*G
        end do
        
        do j = 1, n_y-1
            k = (n_x-1)*(j-1) + 1
            F(k,1) = F(k,1) - A_x*G
        end do
    end subroutine generate_actions_matrix

    subroutine  generate_stiffness_matrix()
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
    end subroutine generate_stiffness_matrix

    subroutine matrices_to_csv()
        101 format(1x, *(g0, :, ", ")) 
    
    ! Open connection (i.e. create file where to write)

        open(unit = 10, access = "sequential", action = "write", &
         status = "replace", file = "S_matrix.csv", form = "formatted") 

    ! Loop across rows

        do k=1, p
            write(10, 101) S(k,:)
        end do 
    ! Close connection
        close(10)

        open(unit = 11, access = "sequential", action = "write", &
         status = "replace", file = "F_matrix.csv", form = "formatted") 

    ! Loop across rows

        do k=1, p
            write(11, 101) F(k,:)
        end do 
    ! Close connection
        close(11)
    end subroutine matrices_to_csv
end subroutine generate_matrices