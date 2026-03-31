program ftcs
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
    
    ! Matrizes de temperatura
    real(8), allocatable :: T(:, :), T_new(:, :), Q_gauss(:, :)

    ! Variáveis auxiliares
    integer :: i, j, n
    real(8) :: rx, ry
    real(8) :: x, y, tempo
    real(8) :: Q, S
    real(8) :: norma

    ! Variáveis para análises
    real(8) :: Tmax_atual
    real(8) :: t_inicio_cpu, t_fim_cpu, tempo_cpu
    logical :: atingiu_limiar = .false.
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

    ! Propriedades Termofísicas e Discretização
    alpha = k / (rho * cp)
    dx = Lx / dble(Nx - 1)
    dy = Ly / dble(Ny - 1)
    
    dt = fator_dt * ((dx**2) / (4.0 * alpha))
    Nt = int(tmax / dt)

    rx = alpha * dt / (dx**2) 
    ry = alpha * dt / (dy**2)

    ! Condição inicial
    T = Tinf 
    T_new = Tinf 
    T(:,1) = Tinf

    ! Pré cálculo do hotspot (Gaussiana)
    do j = 1, Ny
        y = (j-1)*dy
        do i = 1, Nx
            x = (i-1)*dx
            Q_gauss(i,j) = Q0 * exp(-((x-x0)**2 + (y-y0)**2)/(2.0*sigma**2))
        end do
    end do

    ! Abertura dos arquivos de histórico no tempo
    open(20, file = "serie_tempo_ftcs.dat", status="replace")
    open(60, file = "norma_transiente_ftcs.dat", status="replace")

    ! Cronômetro do processador
    call cpu_time(t_inicio_cpu)

    ! Calcula o intervalo de tempo para salvar 4 snapshots equidistantes
    intervalo_snap = tmax / 4.0

    ! ===================================================================
    ! Loop temporal (Cálculos Principais)
    ! ===================================================================
    do n = 1, Nt
        tempo = n * dt
        norma = 0.0 

        ! Duty-Cycle
        if ( mod(tempo, periodo) < ((1.0 - duty) * periodo) ) then
            S = 0.0  ! Começa DESLIGADA
        else
            S = 1.0  ! Termina LIGADA
        end if

        ! Cálculo do FTCS nos nós internos
        do j = 2, Ny-1
            do i = 2, Nx-1
                Q = Q_gauss(i,j) * S

                ! Equação do FTCS
                T_new(i, j) = T(i, j) &
                + rx * (T(i+1, j) - 2.0*T(i, j) + T(i-1, j)) &
                + ry * (T(i, j+1) - 2.0*T(i, j) + T(i, j-1)) &
                + dt * Q / (rho*cp)
                
                ! Variação física máxima de temperatura neste passo
                norma = max(norma, abs(T_new(i, j) - T(i, j)))
            end do
        end do

        ! Condições de Contorno
        do j = 1, Ny
            T_new(1,j)  = T_new(2,j)      ! Esquerda adiabática
            T_new(Nx,j) = T_new(Nx-1,j)   ! Direita adiabática
        end do
        T_new(:,1) = Tinf                 ! Base fixa (Dirichlet)
        do i = 1, Nx
            T_new(i,Ny) = (k/dy*T_new(i,Ny-1) + h*Tinf)/(k/dy + h) ! Topo convectivo (Robin)
        end do

        ! Avança no tempo
        T = T_new
        write(60, *) tempo, norma

        Tmax_atual = maxval(T)
        
        ! Trava de segurança para instabilidade
        if (Tmax_atual > 500.0) exit
        
        ! Verificação do limiar térmico
        if (.not. atingiu_limiar .and. Tmax_atual >= temp_alvo) then
            atingiu_limiar = .true.
            tempo_limiar = tempo
        end if
        
        write(20, *) tempo, T(Nx/2, Ny/2), Tmax_atual, norma

        ! Snapshots adaptativos
        if (tempo >= (proximo_snap * intervalo_snap) .and. proximo_snap <= 4) then
            write(nome_arquivo, '("snap_ftcs_fase", I2.2, ".dat")') proximo_snap
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

    ! Para o cronômetro
    call cpu_time(t_fim_cpu)
    tempo_cpu = t_fim_cpu - t_inicio_cpu

    ! -------------------------------------------------------------------
    ! SALVAMENTO DOS RESULTADOS
    ! -------------------------------------------------------------------
    ! Dados do Mapa de Calor 2D Final
    open(10, file = "ftcs.dat", status="replace")
    do i = 1, Nx
        do j = 1, Ny
            write(10,*) (i-1)*dx, (j-1)*dy, T(i,j)
        end do
    end do
    close(10)
    close(20)

    ! Convergência Espacial
    open(40, file = "dados_convergencia.dat", position="append", status="unknown")
    write(40, *) Nx, dx, maxval(T)
    close(40)

    ! LOG Estruturado
    open(50, file = "log_ftcs.txt", status="replace")
    write(50,'(A)')         "======================================================="
    write(50,'(A)')         "         SIMULADOR TERMICO 2D - FTCS"
    write(50,'(A)')         "======================================================="
    write(50,'(A, I4, A, I4)')    " Malha            : ", Nx, " x ", Ny
    write(50,'(A, F10.6, A)')     " dx               : ", dx, " m"
    write(50,'(A, F10.6, A)')     " dt               : ", dt, " s"
    write(50,'(A, I8)')           " Nt               : ", Nt
    write(50,'(A, F10.3)')        " rx               : ", rx
    write(50,'(A, F10.3, A)')     " CPU              : ", tempo_cpu, " s"
    write(50,'(A, F10.2, A)')     " T_max final      : ", maxval(T), " °C"
    write(50,'(A, F10.3)')        " Media Iter/Passo : ", 1.0d0
    write(50,'(A, ES10.2)')       " Residuo Maximo   : ", 0.0d0
    write(50,'(A, ES10.2)')       " Norma            : ", 0.0d0
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

    ! Dados brutos para tabela (MATLAB)
    open(51, file = "metricas_ftcs.dat", status="replace")
    if (atingiu_limiar) then
        write(51, *) Nx, dx, dt, rx, tempo_cpu, maxval(T), tempo_limiar, 1.0d0, 0.0d0, 0.0d0, Q0, x0, y0, duty
    else
        write(51, *) Nx, dx, dt, rx, tempo_cpu, maxval(T), -1.0d0, 1.0d0, 0.0d0, 0.0d0, Q0, x0, y0, duty
    end if
    close(51)

    deallocate( T, T_new, Q_gauss )

end program ftcs