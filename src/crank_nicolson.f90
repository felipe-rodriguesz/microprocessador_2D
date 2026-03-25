program crank_nicolson
    implicit none

    ! -------------------------------------------------------------------
    ! DECLARAÇÕES DE VARIÁVEIS E CONSTANTES
    ! -------------------------------------------------------------------

    ! Propriedades físicas e Convectivas
    real(8), parameter :: Lx = 0.02, Ly = 0.02
    real(8), parameter :: rho = 2330.0, cp = 700.0, k = 150.0
    real(8), parameter :: alpha = k/(rho*cp)
    real(8), parameter :: h = 10.0, Tinf = 27.0
    
    ! Hotspot
    real(8) :: Q0, x0, y0
    real(8), parameter :: sigma = 0.002, periodo = 0.5
    real(8) :: duty
    real(8) :: S_next, Q_next

    ! Parâmetros da malha e tempo
    integer :: Nx, Ny
    real(8) :: dx, dy, dt, tmax, fator_dt
    integer :: Nt
    
    ! Matrizes
    real(8), allocatable :: T(:, :), T_new(:, :), Q_gauss(:, :), RHS(:, :)

    ! Variáveis auxiliares
    integer :: i, j, n, iter, total_iter = 0
    real(8) :: rx, ry, x, y, tempo, Q, S, Fx, Fy

    ! Variáveis Gauss-Seidel
    integer, parameter :: max_iter = 5000
    real(8) :: T_old_iter, tol, den
    real(8) :: residuo_max, residuo_local, vizinhos

    real(8) :: norma

    ! Variáveis para análises
    real(8) :: Tmax_atual, t_inicio_cpu, t_fim_cpu, tempo_cpu
    logical :: atingiu_limiar = .false.
    real(8) :: tempo_limiar
    real(8), parameter :: temp_alvo = 29.0

    ! Variáveis para os Snapshots
    integer :: proximo_segundo = 1
    character(len=30) :: nome_arquivo

    ! -------------------------------------------------------------------
    ! ÁREA DE EXECUÇÃO (Cálculos e Loops)
    ! -------------------------------------------------------------------
    write(*,*) "======================================================="
    write(*,*) " Entrada: Nx Ny fator_dt Q0 x0 y0 duty"
    write(*,*) " Padrao : 100 100 1.0 1e8 0.01 0.01 0.5"
    write(*,*) "======================================================="
    read(*,*) Nx, Ny, fator_dt, Q0, x0, y0, duty
    allocate( T(Nx, Ny), T_new(Nx, Ny), Q_gauss(Nx, Ny), RHS(Nx, Ny) )

    dx = Lx / dble(Nx - 1)
    dy = Ly / dble(Ny - 1)
    
    dt = fator_dt * ((dx**2) / (4.0 * alpha))
    tmax = 4.75
    Nt = int(tmax / dt)

    rx = alpha * dt / (dx**2) 
    ry = alpha * dt / (dy**2)

    Fx  = rx / 2.0
    Fy  = ry / 2.0
    tol = 1.0e-6 
    den = 1.0 + 2.0*Fx + 2.0*Fy

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

    ! Impressões de Setup
    write(*,'(A)') "======================================================="
    write(*,'(A)') "        SIMULADOR TERMICO 2D - CRANK-NICOLSON          "
    write(*,'(A)') "======================================================="
    write(*,'(A)') "[PARAMETROS DA MALHA]"
    write(*,'(A, I4, A, I4)')       " Malha Espacial (Nx x Ny) : ", Nx, " x ", Ny
    write(*,'(A, F10.6, A)')        " Passo Espacial (dx)      : ", dx, " m"
    write(*,'(A, F10.6, A)')        " Passo de Tempo (dt)      : ", dt, " s"
    write(*,'(A, I8)')              " Numero de Passos (Nt)    : ", Nt
    write(*,'(A)') ""
    write(*,'(A)') "[PROPRIEDADES FISICAS E NUMERICAS]"
    write(*,'(A, ES10.2, A)')       " Difusividade (alpha)     : ", alpha, " m^2/s"
    write(*,'(A, F10.3)')           " Fator de Estabilidade rx : ", rx
    write(*,'(A)') "[PARAMETROS DE SENSIBILIDADE]"
    write(*,'(A, ES10.2, A)')       " Intensidade Q0           : ", Q0, " W/m^3"
    write(*,'(A, F10.5, A)')        " Centro Hotspot x0        : ", x0, " m"
    write(*,'(A, F10.5, A)')        " Centro Hotspot y0        : ", y0, " m"
    write(*,'(A, F10.3)')           " Duty-Cycle               : ", duty
    write(*,'(A)') "======================================================="
    write(*,'(A)') " Iniciando loop temporal... aguarde."
    write(*,'(A)') ""

    ! Salva o histórico no tempo
    ! col 1: tempo | col 2: T_centro | col 3: T_max | col 4: norma
    open(20, file = "serie_tempo_cn.dat", status="replace")

    ! Norma por iteração do passo n=1 (para gráfico de convergência)
    open(60, file = "residuo_iter_cn.dat", status="replace")

    call cpu_time(t_inicio_cpu)

    ! ===================================================================
    ! Loop temporal
    ! ===================================================================
    do n = 1, Nt
        tempo = n * dt

        ! Duty-Cycle no instante t
        if ( mod(tempo, periodo) < (duty * periodo) ) then
            S = 1.0
        else
            S = 0.0
        end if

        ! Duty-Cycle no instante t+dt
        if ( mod(tempo + dt, periodo) < (duty * periodo) ) then
            S_next = 1.0
        else
            S_next = 0.0
        end if

        ! Lado direito — parte explícita do Crank-Nicolson
        do j = 2, Ny-1
            do i = 2, Nx-1
                Q      = Q_gauss(i,j) * S
                Q_next = Q_gauss(i,j) * S_next
                
                RHS(i,j) = T(i,j) &
                         + Fx * (T(i+1,j) - 2.0*T(i,j) + T(i-1,j)) &
                         + Fy * (T(i,j+1) - 2.0*T(i,j) + T(i,j-1)) &
                         + (dt * (Q + Q_next)) / (2.0 * rho * cp)
            end do
        end do

        T_new = T                              ! Chute inicial

        ! --- SOLVER ITERATIVO DE GAUSS-SEIDEL ---
        do iter = 1, max_iter
            norma       = 0.0
            residuo_max = 0.0

            do j = 2, Ny-1
                do i = 2, Nx-1

                    T_old_iter = T_new(i,j)

                    vizinhos = Fx * (T_new(i+1, j) + T_new(i-1, j)) &
                             + Fy * (T_new(i, j+1) + T_new(i, j-1))

                    residuo_local = RHS(i,j) + vizinhos - (den * T_old_iter)

                    T_new(i, j) = T_old_iter + (residuo_local / den)

                    ! Norma do professor: max|T_new - T_old|
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

            ! Salva norma por iteração apenas no passo n=1
            if (n == 1) write(60, *) iter, norma

            ! Critério de parada: norma do professor
            if (norma < tol) exit

        end do

        if (iter >= max_iter) print *, "AVISO: Não convergiu no tempo ", tempo
        total_iter = total_iter + iter

        T = T_new

        Tmax_atual = maxval(T)
        
        if (.not. atingiu_limiar .and. Tmax_atual >= temp_alvo) then
            atingiu_limiar = .true.
            tempo_limiar = tempo
        end if

        ! col 4 = norma ao final do solver (já convergida)
        write(20, *) tempo, T(Nx/2, Ny/2), Tmax_atual, norma

        if (tempo >= (proximo_segundo) .and. proximo_segundo <= 4) then
            write(nome_arquivo, '("snap_cn_", I1, "s.dat")') proximo_segundo
            open(30, file=trim(nome_arquivo), status="replace")
            do i = 1, Nx
                do j = 1, Ny
                    write(30,*) (i-1)*dx, (j-1)*dy, T(i,j)
                end do
            end do
            close(30)
            print *, "--> Snapshot termico salvo: ", trim(nome_arquivo)
            proximo_segundo = proximo_segundo + 1
        end if

    end do

    close(60)
    call cpu_time(t_fim_cpu)
    tempo_cpu = t_fim_cpu - t_inicio_cpu

    write(*,'(A)') ""
    write(*,'(A)') "======================================================="
    write(*,'(A)') "                 RELATORIO DE EXECUCAO                 "
    write(*,'(A)') "======================================================="
    write(*,'(A, F10.3, A)')        " Tempo de CPU Consumido   : ", tempo_cpu, " segundos"
    if (atingiu_limiar) then
        write(*,'(A, F6.1, A, F10.3, A)') " Limiar Termico (", temp_alvo, "°C)   : Atingido em ", tempo_limiar, " s"
    else
        write(*,'(A, F6.1, A)')           " Limiar Termico (", temp_alvo, "°C)   : Nao atingido."
    end if
    write(*,'(A, F10.3)')           " Media de Iteracoes/Passo : ", real(total_iter)/real(Nt)
    write(*,'(A, F10.2, A)')        " Temperatura Max Final    : ", maxval(T), " °C"
    write(*,'(A, ES10.2)')          " Residuo Maximo Final     : ", residuo_max
    write(*,'(A, ES10.2)')          " Norma Final (max|dT|)    : ", norma
    write(*,'(A)') "======================================================="

    open(10, file = "crank_nicolson.dat", status="replace")
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

    open(50, file = "log_cn.txt", status="replace")
    write(50,'(A)')         "======================================================="
    write(50,'(A)')         "         SIMULADOR TERMICO 2D - CRANK-NICOLSON"
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

    ! MÉTRICAS PARA MATLAB (14 colunas — col 10 agora é norma max|dT|)
    open(51, file = "metricas_cn.dat", status="replace")
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

    deallocate( T, T_new, Q_gauss, RHS )

end program crank_nicolson