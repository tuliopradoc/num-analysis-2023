! EULER'S METHOD (SYSTEM OF ODE's):
! ___________________________________________________________________________________

! Approximating the soluution of the system of ODE's 

! u' = -4u + 3v + 6         ;     u(0)=0
! v' = -2.4u + 1.6v + 3.6   ;    v(0)=0

! at a finite number of equally spaced points in [0,0.5].
! ___________________________________________________________________________________

program eulerv
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

    integer :: i, k, m, n
    real(dp) :: u, v, a, h, x
    real(dp) :: Y_arr(2,1), C_arr(2,2), D_arr(2,1), matriz_const1(4), matriz_const2(2)
    real(dp), allocatable :: E(:,:), app_arr(:)
    
    
    u = 0; v = 0; a = 0; x = a; h = 0.1d0; m = size(Y_arr); n = 5


    matriz_const1 = [-4d0, -2.4d0, 3d0, 1.6d0]
    matriz_const2 = [6d0, 3.6d0]

    Y_arr = reshape([u,v], shape(Y_arr))
    C_arr = reshape(matriz_const1, shape(C_arr))
    D_arr = reshape(matriz_const2, shape(D_arr))
    

    allocate(E(size(Y_arr, 1), size(C_arr, 2)))
    allocate(app_arr(m*n+1))
    
    app_arr=[0]

    do i = 1, m 
        E = matmul(C_arr, Y_arr) + D_arr
        do k = 1, n 
            Y_arr(i,1) = Y_arr(i,1) + h*E(i,1)
            x = a + real(k, dp)*h
            app_arr = [app_arr, Y_arr(i,1)]
            E(i,:) = matmul(C_arr(i, :), Y_arr) + D_arr(i,:)
        end do  
        Y_arr(i,1) = u
    end do
    
    print *, app_arr
    
end program eulerv