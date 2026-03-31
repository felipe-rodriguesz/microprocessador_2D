program btcs
    implicit none

    ! -------------------------------------------------------------------
    ! DECLARAÇÕES DE VARIÁVEIS E CONSTANTES
    ! -------------------------------------------------------------------
    ! Variáveis Físicas
    real(8) :: Lx, Ly
    real(8) :: rho, cp, k, alpha
    real(8) :: h, Tinf
    real(8) :: sigma, periodo, tmax
    
    ! Hotspot e Malha
    real(8) :: Q0, x0, y0, duty
    integer :: Nx, Ny
    real(8) :: dx, dy, dt, fator_dt
    integer :: Nt
    
    ! Matrizes
    real(8), allocatable :: T(:, :), T_new(:, :), Q_gauss(:, :)

    ! Variáveis auxiliares
    integer :: i, j, n, iter, total_iter = 0
    real(8) :: rx, ry, x, y, tempo, Q, S

    ! Variáveis Gauss-Seidel
    integer, parameter :: max_iter = 5000
    real(8) :: T_old_iter, tol, den
    real(8) :: residuo_max, vizinhos, RHS_local, residuo_local
    real(8) :: norma

    ! Variáveis para análises
    real(8) :: Tmax_atual, t_inicio_cpu, t_fim_cpu, tempo_cpu
    logical :: atingiu_limiar = .false.
    logical :: gravou_iter = .false.
    real(8) :: tempo_limiar
    real(8), parameter :: temp_alvo = 29.0

    ! Variáveis para os Snapshots
    integer :: proximo_snap = 1
    real(8) :: intervalo_snap
    character(len=30) :: nome_arquivo

    ! ===================================================================
    ! LEITURA DO ARQUIVO DE CONFIGURAÇÃO FÍSICA
    ! ===================================================================
    open(15, file="parametros.txt", status="old", action="read")
    read(15, *) Lx, Ly
    read(15, *) rho, cp, k
    read(15, *) h, Tinf
    read(15, *) sigma, periodo, tmax
    read(15, *) Nx, Ny, fator_dt
    read(15, *) Q0, x0, y0, duty
    close(15)

    ! Alocação das matrizes
    allocate( T(Nx, Ny), T_new(Nx, Ny), Q_gauss(Nx, Ny) )

    ! Difusividade
    alpha = k / (rho * cp)
    dx = Lx / dble(Nx - 1)
    dy = Ly / dble(Ny - 1)
    
    dt = fator_dt * ((dx**2) / (4.0 * alpha))
    Nt = int(tmax / dt)

    rx = alpha * dt / (dx**2) 
    ry = alpha * dt / (dy**2) 

    tol = 1.0e-6 
    den = 1.0 + 2.0*rx + 2.0*ry

    ! Condição inicial
    T = Tinf 
    T_new = Tinf 
    T(:,1) = Tinf

    ! Pré cálculo do hotspot
    do j = 1, Ny
        y = (j-1)*dy
        do i = 1, Nx
            x = (i-1)*dx
            Q_gauss(i,j) = Q0 * exp(-((x-x0)**2 + (y-y0)**2)/(2.0*sigma**2))
        end do
    end do

    ! Salva o histórico no tempo
    open(20, file = "serie_tempo_btcs.dat", status="replace")

    ! Norma por iteração (para gráfico de convergência)
    open(60, file = "residuo_iter_btcs.dat", status="replace")

    ! Cronômetro do processador
    call cpu_time(t_inicio_cpu)

    ! Calcula o intervalo de tempo para salvar 4 snapshots equidistantes
    intervalo_snap = tmax / 4.0

    ! ===================================================================
    ! Loop temporal
    ! ===================================================================
    do n = 1, Nt
        tempo = n * dt

        ! Duty-Cycle
        if ( mod(tempo, periodo) < ((1.0 - duty) * periodo) ) then
            S = 0.0  ! Começa DESLIGADA
        else
            S = 1.0  ! Termina LIGADA
        end if

        T_new = T

        ! --- SOLVER ITERATIVO DE GAUSS-SEIDEL ---
        do iter = 1, max_iter
            norma = 0.0
            residuo_max = 0.0

            do j = 2, Ny-1
                do i = 2, Nx-1

                    Q = Q_gauss(i,j) * S
                    T_old_iter = T_new(i,j)

                    vizinhos  = rx * (T_new(i+1, j) + T_new(i-1, j)) &
                              + ry * (T_new(i, j+1) + T_new(i, j-1))

                    RHS_local = T(i,j) + (dt * Q) / (rho*cp)
                    residuo_local = RHS_local + vizinhos - (den * T_old_iter)

                    T_new(i, j) = T_old_iter + (residuo_local / den)

                    norma = max(norma, abs(T_new(i,j) - T_old_iter))

                    if (abs(residuo_local) > residuo_max) residuo_max = abs(residuo_local)

                end do
            end do

            ! Condições de Contorno
            do j = 1, Ny
                T_new(1,j)  = T_new(2,j)
                T_new(Nx,j) = T_new(Nx-1,j)
            end do
            T_new(:,1) = Tinf
            do i = 1, Nx
                T_new(i,Ny) = (k/dy*T_new(i,Ny-1) + h*Tinf)/(k/dy + h)
            end do

            ! Salva norma no primeiro passo em que a fonte LIGA
            if (S == 1.0 .and. .not. gravou_iter) write(60, *) iter, norma

            ! Critério de parada
            if (norma < tol) exit

        end do

        ! Trava o gravador para não salvar ficheiros gigantes
        if (S == 1.0) gravou_iter = .true.

        total_iter = total_iter + iter
        T = T_new
        Tmax_atual = maxval(T)
        
        if (.not. atingiu_limiar .and. Tmax_atual >= temp_alvo) then
            atingiu_limiar = .true.
            tempo_limiar = tempo
        end if

        write(20, *) tempo, T(Nx/2, Ny/2), Tmax_atual, norma

        ! Snapshots adaptativos
        if (tempo >= (proximo_snap * intervalo_snap) .and. proximo_snap <= 4) then
            write(nome_arquivo, '("snap_btcs_fase", I2.2, ".dat")') proximo_snap
            open(30, file=trim(nome_arquivo), status="replace")
            do i = 1, Nx
                do j = 1, Ny
                    write(30,*) (i-1)*dx, (j-1)*dy, T(i,j)
                end do
            end do
            close(30)
            proximo_snap = proximo_snap + 1
        end if

    end do

    close(60)
    call cpu_time(t_fim_cpu)
    tempo_cpu = t_fim_cpu - t_inicio_cpu

    ! -------------------------------------------------------------------
    ! SALVAMENTO DOS RESULTADOS E MÉTRICAS
    ! -------------------------------------------------------------------
    open(10, file = "backward_euler.dat", status="replace")
    do i = 1, Nx
        do j = 1, Ny
            write(10,*) (i-1)*dx, (j-1)*dy, T(i,j)
        end do
    end do
    close(10)
    close(20)

    open(40, file = "dados_convergencia.dat", position="append", status="unknown")
    write(40, *) Nx, dx, maxval(T)
    close(40)

    open(50, file = "log_btcs.txt", status="replace")
    write(50,'(A)')         "======================================================="
    write(50,'(A)')         "         SIMULADOR TERMICO 2D - BACKWARD EULER"
    write(50,'(A)')         "======================================================="
    write(50,'(A, I4, A, I4)')    " Malha            : ", Nx, " x ", Ny
    write(50,'(A, F10.6, A)')     " dx               : ", dx, " m"
    write(50,'(A, F10.6, A)')     " dt               : ", dt, " s"
    write(50,'(A, I8)')           " Nt               : ", Nt
    write(50,'(A, F10.3)')        " rx               : ", rx
    write(50,'(A, F10.3, A)')     " CPU              : ", tempo_cpu, " s"
    write(50,'(A, F10.2, A)')     " T_max final      : ", maxval(T), " °C"
    write(50,'(A, F10.3)')        " Media Iter/Passo : ", real(total_iter)/real(Nt)
    write(50,'(A, ES10.2)')       " Residuo Maximo   : ", residuo_max
    write(50,'(A, ES10.2)')       " Norma (max|dT|)  : ", norma
    if (atingiu_limiar) then
        write(50,'(A, F10.3, A)') " Limiar (29°C)    : ", tempo_limiar, " s"
    else
        write(50,'(A)')           " Limiar (29°C)    : Nao atingido"
    end if
    write(50,'(A)')         "[PARAMETROS DE SENSIBILIDADE]"
    write(50,'(A, ES10.2, A)')    " Q0               : ", Q0, " W/m^3"
    write(50,'(A, F10.5, A)')     " x0               : ", x0, " m"
    write(50,'(A, F10.5, A)')     " y0               : ", y0, " m"
    write(50,'(A, F10.3)')        " Duty-Cycle       : ", duty
    write(50,'(A)')         "======================================================="
    close(50)

    ! Métricas para o MATLAB
    open(51, file = "metricas_btcs.dat", status="replace")
    if (atingiu_limiar) then
        write(51, *) Nx, dx, dt, rx, tempo_cpu, maxval(T), tempo_limiar, &
                     real(total_iter)/real(Nt), residuo_max, norma, &
                     Q0, x0, y0, duty
    else
        write(51, *) Nx, dx, dt, rx, tempo_cpu, maxval(T), -1.0d0, &
                     real(total_iter)/real(Nt), residuo_max, norma, &
                     Q0, x0, y0, duty
    end if
    close(51)

    deallocate( T, T_new, Q_gauss )

end program btcs