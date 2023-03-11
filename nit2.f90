program nit
    implicit none

    integer :: n, i, j, k, l
    real(8) :: N_cells, I_cells, T_cells, h, x, g, p
    real(8) :: NIT_arr(3), NIT_matrix(3,1), NIT_matrix_transp(1,3), c1x3(1,3), c1x1(1,1), c3x3(3,3), temp1(3,1), NIT_temp(3,1)
    real(8) :: NIT_3x3(3,3,3), NIT_1x3(3,1,3), NIT_1x1(3,1,1), f1(3,1,1), f2(3,1,1), temp2(3,1), M(3,1)
    real(8) :: N_arr_3x3(9), I_arr_3x3(9), T_arr_3x3(9), N_arr_1x3(3), I_arr_1x3(3), T_arr_1x3(3) 
    real(8), allocatable :: NIT_matrices(:), NIT_csv(:,:)
    
    N_cells = 1d0; I_cells = 1.22d0; T_cells = 1d0
    
    h = 0.125d0        ! step
    n = 240            ! iterations
    
    NIT_arr = [N_cells, I_cells, T_cells]
    
    allocate(NIT_csv(3,n + 1))
    allocate(NIT_matrices(3*n + 3))
    NIT_matrices = NIT_arr                              ! array to create .csv file

    NIT_matrix = reshape(NIT_arr, shape(NIT_matrix))                    ! 3x1 matrix
    NIT_matrix_transp = transpose(NIT_matrix)                           ! 1x3 matrix

    g = (0.2710d0)/(1.8130d0) - 0.8130d0
    p = (0.7829d0)/(1.8620d0) - 0.3634d0

    ! write all columns in a row to apply reshape

    N_arr_3x3 = [real(-1.289288e-6, kind = 8), -0.1379d0, -0.9314d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0]      
    I_arr_3x3 = [0d0, g, 0d0, 0d0, 0d0, p, 0d0, 0d0, 0d0]
    T_arr_3x3 = [0d0, 0d0, 1.1890d0, 0d0, 0d0, -0.1469d0, 0d0, 0d0, -0.17704d0]

    N_arr_1x3 = [0.4312d0, 0d0, 0d0]
    I_arr_1x3 = [0d0, -0.57d0, 0d0]
    T_arr_1x3 = [0d0, 0d0, 0.4426d0]
    
    NIT_3x3(1,:,:) = reshape(N_arr_3x3, shape(c3x3))
    NIT_3x3(2,:,:) = reshape(I_arr_3x3, shape(c3x3))        ! 3x3 coefficient matrices
    NIT_3x3(3,:,:) = reshape(T_arr_3x3, shape(c3x3))

    NIT_1x3(1,:,:) = reshape(N_arr_1x3, shape(c1x3))
    NIT_1x3(2,:,:) = reshape(I_arr_1x3, shape(c1x3))        ! 1x3  coefficient matrices
    NIT_1x3(3,:,:) = reshape(T_arr_1x3, shape(c1x3))
    
    NIT_1x1(1,:,:) = reshape([0d0], shape(c1x1))
    NIT_1x1(2,:,:) = reshape([0.7d0], shape(c1x1))          ! 1x1 constant matrices 
    NIT_1x1(3,:,:) = reshape([0d0], shape(c1x1))

    ! Applying Runge-Kutta method of order 2

    do i = 1, n
        do j = 1, 3
            temp1 = matmul(NIT_3x3(j,:,:), NIT_matrix)
            f1(j,:,:) = matmul(NIT_matrix_transp, temp1) + matmul(NIT_1x3(j,:,:), NIT_matrix) + NIT_1x1(j,:,:)
            M(j,1) = NIT_matrix(j,1) + (h/2d0)*f1(j,1,1)
        end do
        do l = 1, 3
            temp2 = matmul(NIT_3x3(l,:,:), M)
            f2(l,:,:) = matmul(NIT_matrix_transp, temp2) + matmul(NIT_1x3(l,:,:), M) + NIT_1x1(l,:,:)
            NIT_temp(l,1) = NIT_matrix(l,1) + h*f2(l,1,1)
            NIT_matrices = [NIT_matrices, NIT_temp(l,1)]
        end do
        NIT_matrix = NIT_temp
        NIT_matrix_transp = transpose(NIT_temp)
        NIT_3x3(2,2,1) = 0.2710d0/(0.8130d0+NIT_matrix(1,1)) - 0.8130d0            ! calculating p and g to use in the next iteration
        NIT_3x3(2,3,2) = 0.7829d0/(0.8620d0+NIT_matrix(3,1)) - 0.3634d0
        x = h*real(i, kind = 8)
    end do
    
    NIT_csv = reshape(NIT_matrices, shape(NIT_csv))         ! matrix to generate the .csv file

    101 format(1x, *(g0, :, ", ")) 
    
    ! Open connection (i.e. create file where to write)

    open(unit = 10, access = "sequential", action = "write", &
         status = "replace", file = "NIT_RK.csv", form = "formatted") 

    ! Loop across rows

    do k=1,3
       write(10, 101) NIT_csv(k,:)
    end do 
    ! Close connection
    close(10)

    deallocate(NIT_matrices)
    deallocate(NIT_csv)
end program nit