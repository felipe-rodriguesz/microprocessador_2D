! =============================================================
! Simulador térmico 2D - Euler Implícito + Gauss-Seidel
! =============================================================

program simulador_chip
    implicit none

    ! Parâmetros de entrada e derivados
    real(8) :: Lx, Ly, k, rho, cp, alpha, Q0, x0, y0, sigma, T0, Tinf, h
    real(8) :: dt, tmax, t_ciclo_ini, t_ciclo_fim, period, duty, t_on, tol, dx, dy, aE, aW, aN, aS, aP
    integer :: Nx, Ny, Nt, max_iter
    character(len=*), parameter :: out_dir = '../data/results/'

    real(8) :: tempo, St, T_new, R, R_sum, norm_R, norm_R0, residuo_norm, N_pontos, T_max_atual
    integer :: i, j, n, iter
    integer :: iter_conv
    real(8), allocatable :: T(:,:), Told(:,:), b(:,:), Q(:,:)

    call init_parametros()
    allocate(T(Nx, Ny), Told(Nx, Ny), b(Nx, Ny), Q(Nx, Ny))
    
    N_pontos = dble((Nx - 2) * (Ny - 2))
    T = T0 ! Chute inicial uniforme

    ! Termo fonte
    do j = 1, Ny
        do i = 1, Nx
            Q(i,j) = Q0 * exp(-(((i-1)*dx-x0)**2 + ((j-1)*dy-y0)**2) / (2.0d0*sigma**2))
        end do
    end do

    open(10, file=trim(out_dir)//'norma_iter_chip.dat', status='replace')
    open(11, file=trim(out_dir)//'convergencia_passos.dat', status='replace')
    open(12, file=trim(out_dir)//'tmax_tempo.dat',          status='replace')
    open(13, file=trim(out_dir)//'sensores_tempo.dat',      status='replace')

    ! Loop Temporal
    do n = 1, Nt
        tempo = dble(n) * dt
        Told = T
        
        ! Função temporal (St)
        St = 0.0d0
        if (tempo >= t_ciclo_ini .and. tempo <= t_ciclo_fim) then
            if (mod(tempo - t_ciclo_ini, period) <= t_on) St = 1.0d0
        end if

        ! Termo fonte
        do j = 2, Ny-1
            do i = 2, Nx-1
                b(i,j) = Told(i,j)/dt + (Q(i,j)/(rho*cp)) * St
            end do
        end do

        T(:,1)  = T0                                     ! Base (Dirichlet)

        ! Gauss-Seidel
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

            ! Condições de Contorno
            T(1,:)  = T(2,:)                                    ! Parede Esquerda Isolada (Neumann)
            T(Nx,:) = T(Nx-1,:)                                 ! Parede Direita Isolada (Neumann)
            do i = 1, Nx
                T(i,Ny) = (k*T(i,Ny-1)/dy + h*Tinf)/(k/dy+h)    ! Topo com Convecção (Robin)
            end do

            ! Resíduo
            norm_R = sqrt(R_sum / N_pontos)
            if (iter == 1) norm_R0 = max(norm_R, 1.0d-16)
            residuo_norm = norm_R / norm_R0

            if (n == 1) write(10,*) iter, residuo_norm
            if (residuo_norm < tol) exit
        end do
        iter_conv = iter   ! guarda quantas iterações foram necessárias

        write(11,*) n, tempo, iter_conv, residuo_norm

        T_max_atual = maxval(T)
        write(12,*) n, tempo, T_max_atual, St
        write(13,*) tempo, T(Nx/2, Ny/2), T(Nx/4, Ny/2), T(1, Ny/2)
        write(*,'(A,I5,A,I5,A,ES10.3)') "Passo: ", n, " | Iteracoes: ", iter, " | Residuo: ", residuo_norm
    end do
    close(10)
    close(11)
    close(12)
    close(13)

    ! Mapa de Calor
    open(50, file=trim(out_dir)//'mapa_calor_chip.dat', status='replace')
    do j = 1, Ny
        do i = 1, Nx
            write(50, '(2F12.6, ES15.6)') (i-1)*dx, (j-1)*dy, T(i,j)
        end do
        write(50, *)
    end do
    close(50)

    print *, "Simulacao do Microprocessador Concluida!"
    print *, "O Processador atingiu: ", maxval(T), " C"
    deallocate(T, Told, b, Q)

contains

    subroutine init_parametros()
        namelist /parametros/ Lx, Ly, Nx, Ny, k, rho, cp, Q0, x0, y0, sigma, T0, Tinf, h, dt, tmax, &
                      t_ciclo_ini, t_ciclo_fim, period, duty, tol, max_iter

        open(99, file='../parametros.txt', status='old')
        read(99, nml=parametros)
        close(99)

        ! Cálculos derivados e coeficientes
        alpha = k / (rho * cp)
        t_on  = duty * period
        Nt    = int(tmax / dt)
        dx    = Lx / dble(Nx - 1)
        dy    = Ly / dble(Ny - 1)

        aE = alpha / dx**2
        aW = alpha / dx**2
        aN = alpha / dy**2
        aS = alpha / dy**2
        aP = 1.0d0/dt + aE + aW + aN + aS

        call execute_command_line('mkdir -p ' // trim(out_dir))
    end subroutine init_parametros

end program simulador_chip