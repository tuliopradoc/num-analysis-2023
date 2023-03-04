program nit
    implicit none

    integer :: n, i, j
    real(8) :: N_cells, I_cells, T_cells, h, x, g, p
    real(8) :: NIT_arr(3), NIT_matrix(3,1), NIT_matrix_transp(1,3), c1x3(1,3), c1x1(1,1), c3x3(3,3), temp(3,1), NIT_temp(3,1)
    real(8) :: NIT_3x3(3,3,3), NIT_1x3(3,1,3), NIT_1x1(3,1,1), f(3,1,1), NIT_matrices(240,3,1)
    real(8) :: N_arr_3x3(9), I_arr_3x3(9), T_arr_3x3(9), N_arr_1x3(3), I_arr_1x3(3), T_arr_1x3(3) 
    
    N_cells = 1d0; I_cells = 1.22d0; T_cells = 1d0
    
    h = 0.125d0        ! step
    n = 240            ! iterations
    
    NIT_arr = [N_cells, I_cells, T_cells]

    NIT_matrix = reshape(NIT_arr, shape(NIT_matrix))                    ! 3x1 matrix
    NIT_matrix_transp = transpose(NIT_matrix)                           ! 1x3 matrix

    g = (0.2710d0)/(1.8130d0) - 0.8130d0
    p = (0.7829d0)/(1.8620d0) - 0.3634d0

    N_arr_3x3 = [real(-1.289288e-6, kind = 8), -0.1379d0, -0.9314d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0]      ! write all columns in a row to apply reshape
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


    do i = 1, n
        do j = 1, 3
            temp = matmul(NIT_3x3(j,:,:), NIT_matrix)
            f(j,:,:) = matmul(NIT_matrix_transp, temp) + matmul(NIT_1x3(j,:,:), NIT_matrix) + NIT_1x1(j,:,:)
            NIT_temp(j,1) = NIT_matrix(j,1) + h*f(j,1,1)
        end do
        NIT_matrix = NIT_temp
        NIT_matrix_transp = transpose(NIT_temp)
        NIT_matrices(i,:,:) = NIT_matrix
        NIT_3x3(2,2,1) = 0.2710d0/(0.8130d0+NIT_matrix(1,1)) - 0.8130d0            ! calculating p and g to use in the next iteration
        NIT_3x3(2,3,2) = 0.7829d0/(0.8620d0+NIT_matrix(3,1)) - 0.3634d0
        x = h*real(i, kind = 8)
    end do

    print *, NIT_matrices(1,:,:)
    print *, NIT_matrices(2,:,:)
    print *, NIT_matrices(240,:,:)
end program nit