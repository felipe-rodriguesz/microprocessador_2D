! =============================================================
! Simulador termico 2D - Regime estacionario + Gauss-Seidel
! =============================================================

program chip_estacionario
    implicit none

    ! Parametros de entrada e derivados
    real(8) :: Lx, Ly, k, rho, cp, alpha, Q0, x0, y0, sigma, T0, Tinf, h
    real(8) :: tol, dx, dy, aE, aW, aN, aS, aP
    integer :: Nx, Ny, max_iter
    character(len=*), parameter :: out_dir = '../data/results/'

    real(8), parameter :: pi = 3.141592653589793d0
    real(8) :: T_new, R, R_sum, norm_R, norm_R0, residuo_norm, N_pontos
    real(8) :: freq_x, freq_y, pos_x, pos_y
    integer :: i, j, iter
    real(8), allocatable :: T(:,:), b(:,:), Q(:,:)

    call init_parametros()
    allocate(T(Nx, Ny), b(Nx, Ny), Q(Nx, Ny))

    N_pontos = dble((Nx - 2) * (Ny - 2))
    ! Chute inicial: senoide de alta frequencia (modo k=30, m=30)
    freq_x = 30.0d0 * pi / Lx
    freq_y = 30.0d0 * pi / Ly
    do j = 1, Ny
        do i = 1, Nx
            pos_x = (i-1)*dx
            pos_y = (j-1)*dy
            T(i,j) = 50.0d0 * sin(freq_x * pos_x) * sin(freq_y * pos_y)
        end do
    end do

    ! Termo fonte espacial
    do j = 1, Ny
        do i = 1, Nx
            Q(i,j) = Q0 * exp(-(((i-1)*dx-x0)**2 + ((j-1)*dy-y0)**2) / (2.0d0*sigma**2))
        end do
    end do

    ! Monta termo independente do sistema estacionario
    do j = 2, Ny-1
        do i = 2, Nx-1
            b(i,j) = Q(i,j)/(rho*cp)
        end do
    end do

    open(10, file=trim(out_dir)//'norma_iter_chip.dat', status='replace')

    ! Solver Gauss-Seidel
    do iter = 1, max_iter
        R_sum = 0.0d0

        do j = 2, Ny-1
            do i = 2, Nx-1
                T_new = (aE*T(i+1,j) + aW*T(i-1,j) + aN*T(i,j+1) + aS*T(i,j-1) + b(i,j)) / aP
                R = T_new - T(i,j)
                R_sum = R_sum + R**2
                T(i,j) = T_new
            end do
        end do

        ! Condicoes de contorno
        T(:,1)  = T0                                     ! Base (Dirichlet)
        T(1,:)  = T(2,:)                                 ! Parede esquerda isolada (Neumann)
        T(Nx,:) = T(Nx-1,:)                              ! Parede direita isolada (Neumann)
        do i = 1, Nx
            T(i,Ny) = (k*T(i,Ny-1)/dy + h*Tinf)/(k/dy+h) ! Topo com conveccao (Robin)
        end do

        norm_R = sqrt(R_sum / N_pontos)
        if (iter == 1) norm_R0 = max(norm_R, 1.0d-16)
        residuo_norm = norm_R / norm_R0
        write(10,*) iter, residuo_norm

        if (residuo_norm < tol) exit
    end do
    write(*,'(A,I5,A,ES10.3)') "Iteracoes GS: ", iter, " | Residuo Maximo: ", residuo_norm

    close(10)

    ! Salvar mapa de calor
    open(50, file=trim(out_dir)//'mapa_calor_chip.dat', status='replace')
    do j = 1, Ny
        do i = 1, Nx
            write(50, '(2F12.6, ES15.6)') (i-1)*dx, (j-1)*dy, T(i,j)
        end do
    end do
    close(50)

    print *, 'Simulacao estacionaria do Microprocessador concluida!'
    print *, 'Temperatura maxima: ', maxval(T), ' C'

    deallocate(T, b, Q)

contains

    subroutine init_parametros()
        namelist /parametros/ Lx, Ly, Nx, Ny, k, rho, cp, Q0, x0, y0, sigma, T0, Tinf, h, tol, max_iter

        open(99, file='../parametros_estacionario.txt', status='old')
        read(99, nml=parametros)
        close(99)

        alpha = k / (rho * cp)
        dx    = Lx / dble(Nx - 1)
        dy    = Ly / dble(Ny - 1)

        aE = alpha / dx**2
        aW = alpha / dx**2
        aN = alpha / dy**2
        aS = alpha / dy**2
        aP = aE + aW + aN + aS

        call execute_command_line('mkdir -p ' // trim(out_dir))
    end subroutine init_parametros

end program chip_estacionario
