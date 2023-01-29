! TAYLOR'S METHOD OF DEGREE 2: 
! ___________________________________________________________________________________

! Program to approximate the value of a solution to the IVP 

! y' = Pe^{cx} + Qe^{kx}
! y(a)= y_{0}

! at a finite number of equally spaced points in [a,b], where b > a is a real number.
! ___________________________________________________________________________________

program ytaylor2
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    ! Declaring variable types 

    integer :: n, i
    real(dp) :: P, Q, c, k, y_0, h, a, b, x, y
    real(dp), allocatable :: array_x(:), array_y(:)

    ! Reading values for b and n

    print *, 'Type a value for b and an INTEGER value for n '
    read *, b, n

    allocate(array_x(n+1))
    allocate(array_y(n+1))
    
    ! Variables

    P = 1; Q = 1; c = 1; k = 1; a = 1; y_0 = 2
    
    h = (b-a)/n
    x = a
    y = y_0

    ! Approximating the values y(x_{i})

    array_x = [x]
    array_y = [y]

    do i = 1, n
        y = y + h*T2(x)
        x = x + real(i, dp)*h
        array_x = [array_x, x]
        array_y = [array_y, y]
    end do
    
    ! Printing the results

    write(*,*) array_x
    write(*,*) array_y

    deallocate(array_x)
    deallocate(array_y)

    ! Function to calculate approximate values of " T_{2}(x_{i}, y_{i}, h) "

    contains

    real(dp) function T2(point) 
        use, intrinsic :: iso_fortran_env, only : dp => real64
        implicit none

        real(dp), intent(in) :: point
        real(dp) ::  f1, f2, f3, l, df
    
        f1 = P*real(exp(c*(point-h)), dp) + Q*real(exp(k*(point-h)), dp)
        f2 = P*real(exp(c*(point+h)), dp) + Q*real(exp(k*(point+h)), dp)
        l = 2_dp*h
        df = (f2-f1)/l
        
        f3 = P*real(exp(c*point), dp) + Q*real(exp(k*point), dp)

        T2 = f3 + (h/2)*df
    end function T2
end program ytaylor2

