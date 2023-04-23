program project2
    implicit none
    integer :: n_x, n_y
    real(8) :: a, b, c, l, h, T, Q, M, E, D, G

    interface
        subroutine generate_matrices(n_x,n_y,h,l,a,b,c,T,Q,M,E,D,G)
            implicit none
            integer :: i, j, k, p, w, n, r
            integer, intent(in) :: n_x, n_y
            real(8), intent(in) :: a, b, c, h, l, T, Q, M, E, D, G
            real(8) :: h_x, h_y, Hx, Hy, A_x, B_y, C_xy, D_xy
            real(8), allocatable :: S(:,:), x(:), y(:), F(:,:), U(:), partial_x(:), partial_y(:)
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

    call generate_matrices_and_solve(n_x,n_y,h,l,a,b,c,T,Q,M,E,D,G)
    

end program project2

subroutine generate_matrices_and_solve(n_x,n_y,h,l,a,b,c,T,Q,M,E,D,G)
    implicit none
    integer :: i, j, k, p, w, n, r
    integer, intent(in) :: n_x, n_y
    real(8), intent(in) :: a, b, c, h, l, T, Q, M, E, D, G
    real(8) :: h_x, h_y, Hx, Hy, A_x, B_y, C_xy, D_xy
    real(8), allocatable :: S(:,:), x(:), y(:), F(:,:), U(:), partial_x(:), partial_y(:)
    
    h_x = 2d0*l/real(n_x, kind=8)
    h_y = 2d0*h/real(n_y, kind=8)
    Hx = 2*real(h_x, kind=8)
    Hy = 2*real(h_y, kind=8)
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
    call solve_system_for_U()
    call calculate_partial_x_and_y()
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
    
    subroutine solve_system_for_U()
        !real(8), allocatable :: U(:)
        interface
            subroutine solve_linear_system(A, b, x, n)
                implicit none
                integer, intent(in) :: n
                real(8), intent(in) :: A(n,n), b(n)
                real(8), intent(out) :: x(n)
                integer :: ipiv(n), info
            end subroutine solve_linear_system
        end interface
    
    	! solutions vector/matrix
    	
        allocate(U(p))
    	
    	! solving the linear system using solve_linear_system subroutine
    	
        call solve_linear_system(S, F, U, p)
    	
    end subroutine solve_system_for_U

    subroutine calculate_partial_x_and_y()
        
        allocate(partial_x(p), partial_y(p))

        partial_x = 0d0
        partial_y = 0d0
    
    
        do j = 1, n_y-1
            do k = 1, p
                r = k - (n_x-1)*(j-1)
                
                if (r/=1 .and. r/=n_x-1) then
                    partial_x(k) = (U(k-1)-U(k+1))/Hx
                end if
                
                if ( r == 1 ) then
                    partial_x(k) = (G-U(k+1))/Hx
                end if
                
                if ( r == n_x -1 ) then
                    partial_x(k) = (U(k-1)-G)/Hx
                end if
            end do
        end do 
    
        do i = 1, n_x-1
            do j = 1, n_y-1
                if ( j==1 ) then
                    w = n_x -1 +i
                    partial_y(i) = (U(w)-G)/Hy
                end if
                
                if ( 1<j .and. j<n_y-1 ) then
                    w = (n_x-1)*j + i
                    n = (n_x-1)*(j-2) + i
                    k = (n_x-1)*(j-1) + i
                    partial_x(k) = (U(w)-U(n))/Hy
                end if
                
                if ( j==n_y-1 ) then
                    n = (n_x-1)*(n_y-3) + i
                    k = (n_x-1)*(n_y-2) + i
                    partial_x(k) = (G-U(n))/Hy
                end if
            end do
        end do
    end subroutine calculate_partial_x_and_y

    subroutine matrices_to_csv()
        101 format(1x, *(g0, :, ", ")) 

        open(unit = 10, access = "sequential", action = "write", & 
        status = "replace", file = "stiffness_matrix.csv", form = "formatted") 
        do k=1, p 
            write(10, 101) S(k,:)
        end do 
        close(10)

        open(unit = 11, access = "sequential", action = "write", &
        status = "replace", file = "actions_matrix.csv", form = "formatted") 
        do k=1, p
            write(11, 101) F(k,:)
        end do 
        close(11)
        
        open(unit = 12, access = "sequential", action = "write", &
        status = "replace", file = "solutions.csv", form = "formatted") 
        do k=1, p
            write(12, 101) U(k)
        end do 
        close(12)

        open(unit = 13, access = "sequential", action = "write", &
        status = "replace", file = "partial_x.csv", form = "formatted") 
        do k=1, p
            write(13, 101) partial_x(k)
        end do 
        close(13)
    
        open(unit = 14, access = "sequential", action = "write", &
        status = "replace", file = "partial_y.csv", form = "formatted") 
        do k=1, p
            write(14, 101) partial_y(k)
        end do 
        close(14)
    end subroutine matrices_to_csv
end subroutine generate_matrices_and_solve

subroutine solve_linear_system(A, b, x, n)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: A(n,n), b(n)
  real(8), intent(out) :: x(n)
  integer :: ipiv(n), info

  ! Make a copy of the input matrix A and vector b
  
  real(8) :: A_copy(n,n), b_copy(n)
  A_copy = A
  b_copy = b

  ! Compute the LU factorization of A_copy
  
  call dgetrf(n, n, A_copy, n, ipiv, info)

  if (info /= 0) then
    print*, "Error: Matrix is singular"
    stop
  end if

  ! Solve the system using the LU factorization and backward substitution
  
  call dgetrs('N', n, 1, A_copy, n, ipiv, b_copy, n, info)

  if (info /= 0) then
    print*, "Error: Failed to solve system"
    stop
  end if

  ! Copy the solution to the output vector x
  
  x = b_copy
end subroutine solve_linear_system


