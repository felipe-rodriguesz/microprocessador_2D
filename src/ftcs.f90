program simulacao_termica_ftcs
    implicit none

    ! -------------------------------------------------------------------
    ! DECLARAÇÕES DE VARIÁVEIS E CONSTANTES
    ! -------------------------------------------------------------------

    ! Propriedades físicas
    real(8), parameter :: Lx = 0.02                 ! Comprimento (m)
    real(8), parameter :: Ly = 0.02                 ! Altura (m)
    real(8), parameter :: rho = 2330.0              ! Densidade (kg/m^3)
    real(8), parameter :: cp = 700.0                ! Calor específico (J/(Kg*K))
    real(8), parameter :: k = 150.0                 ! Condutividade (W/(m*K))                          
    real(8), parameter :: alpha = k/(rho*cp)        ! Difusividade (m^2/s)
    
    ! Parâmetros convectivos
    real(8), parameter :: h = 10.0                  ! Convecção (W/(m^2*K))
    real(8), parameter :: Tinf = 300.0              ! Temp do fluido externo (K)
    
    ! Hotspot
    real(8), parameter :: Q0 = 1.0e8                ! Intesidade máx (W/m^3)
    real(8), parameter :: x0 = Lx/2.0               ! Centro X
    real(8), parameter :: y0 = Ly/2.0               ! Centro Y
    real(8), parameter :: sigma = 0.002             ! Tamanho do Hotspot (m)
    real(8), parameter :: periodo = 0.5             ! Pulso

    ! Parâmetros da malha
    integer, parameter :: Nx = 100, Ny = 100         ! N° de nos da malha
    
    ! Variáveis de Malha e Tempo
    real(8) :: dx, dy, dt, tmax
    integer :: Nt
    
    ! Matrizes de temperatura
    real(8) :: T(Nx, Ny), T_new(Nx, Ny), Q_gauss(Nx, Ny)

    ! Variáveis auxiliares
    integer :: i, j, n                                ! Índices
    real(8) :: rx, ry                                 ! Coeficientes FTCS
    real(8) :: x, y, tempo                            ! Coordenadas e tempo
    real(8) :: Q, S                                   ! Fonte e ativação

    ! Variáveis para análises
    real(8) :: Tmax_atual                             ! Pico térmico
    real(8) :: t_inicio_cpu, t_fim_cpu, tempo_cpu     ! Medir custo computacional
    logical :: atingiu_limiar = .false.               ! Flag para o alerta de temperatura
    real(8) :: tempo_limiar                           ! Guarda o instante exato do limiar
    real(8), parameter :: temp_alvo = 302.0           ! Limiar de aquecimento para análise

    ! -------------------------------------------------------------------
    ! ÁREA DE EXECUÇÃO (Cálculos e Loops)
    ! -------------------------------------------------------------------
    dx = Lx / (Nx - 1)
    dy = Ly / (Ny - 1)
    
    dt = (dx**2) / (4.0 * alpha)
    tmax = 4.75
    Nt = int(tmax / dt)

    rx = alpha * dt / (dx**2) 
    ry = alpha * dt / (dy**2)

    if (rx + ry > 0.5) then
        print *, "ERRO: Esquema instável! rx + ry = ", rx + ry
        stop
    end if

    ! Condição inicial
    T = Tinf 
    T_new = Tinf 
    T(:,1) = Tinf

    ! Pré cálculo do hotspot
    do j = 1, Ny
        y = (j-1)*dy
        do i = 1, Nx
            x = (i-1)*dx
            Q_gauss(i,j) = Q0 * exp(-((x-x0)**2 + (y-y0)**2)/(2.0d0*sigma**2))
        end do
    end do

    ! Impressões de Setup
    write(*,'(A)') "======================================================="
    write(*,'(A)') "        SIMULADOR TERMICO 2D - FTCS          "
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
    write(*,'(A)') "======================================================="
    write(*,'(A)') " Iniciando loop temporal... aguarde."
    write(*,'(A)') ""

    ! Salva o histórico no tempo
    open(20, file = "serie_tempo_ftcs.dat", status="replace")

    ! Cronômetro do processador
    call cpu_time(t_inicio_cpu)

    ! ===================================================================
    ! Loop temporal
    ! ===================================================================
    do n = 1, Nt
        tempo = n * dt

        ! Duty-Cycle baseado no tempo físico
        if ( mod(tempo, periodo) < (periodo / 2.0) ) then
            S = 1.0
        else
            S = 0.0
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
            end do
        end do

        ! Condições de Contorno
        ! 1. Laterais adiabáticas (Newmann)
        do j = 1, Ny
            T_new(1,j)  = T_new(2,j)
            T_new(Nx,j) = T_new(Nx-1,j)
        end do

        ! 2. Base com temperatura fixa (Dirichlet)
        T_new(:,1) = Tinf

        ! 3. Topo convectivo (Robin)
        do i = 1, Nx
            T_new(i,Ny) = (k/dy*T_new(i,Ny-1) + h*Tinf)/(k/dy + h)
        end do

        ! Avança no tempo
        T = T_new

        ! Extração de dados a cada passo de tempo
        Tmax_atual = maxval(T)
        
        ! Verificação do limiar térmico
        if (.not. atingiu_limiar .and. Tmax_atual >= temp_alvo) then
            atingiu_limiar = .true.
            tempo_limiar = tempo
        end if
        
        write(20, *) tempo, T(Nx/2, Ny/2), Tmax_atual
    end do

    ! Para o cronômetro
    call cpu_time(t_fim_cpu)
    tempo_cpu = t_fim_cpu - t_inicio_cpu

    ! Imprime os resultados
    write(*,'(A)') ""
    write(*,'(A)') "======================================================="
    write(*,'(A)') "                 RELATORIO DE EXECUCAO                 "
    write(*,'(A)') "======================================================="
    write(*,'(A, F10.3, A)')        " Tempo de CPU Consumido   : ", tempo_cpu, " segundos"
    if (atingiu_limiar) then
        write(*,'(A, F6.1, A, F10.3, A)') " Limiar Térmico (", temp_alvo, "K)   : Atingido em ", tempo_limiar, " s"
    else
        write(*,'(A, F6.1, A)')           " Limiar Térmico (", temp_alvo, "K)   : Nao atingido."
    end if
    write(*,'(A, F10.2, A)')        " Temperatura Max Final    : ", maxval(T), " K"
    write(*,'(A)') "======================================================="
    write(*,'(A)') " --> Sucesso: 'ftcs.dat' e 'serie_tempo_ftcs.dat' gerados."

    ! -------------------------------------------------------------------
    ! SALVAMENTO DOS RESULTADOS
    ! -------------------------------------------------------------------
    open(10, file = "ftcs.dat", status="replace")

    do i = 1, Nx
        do j = 1, Ny
            write(10,*) (i-1)*dx, (j-1)*dy, T(i,j)
        end do
    end do

    close(10)
    close(20)

end program simulacao_termica_ftcs