program simulacao_termica_btcs
    implicit none

    ! -------------------------------------------------------------------
    ! DECLARAÇÕES DE VARIÁVEIS E CONSTANTES
    ! -------------------------------------------------------------------

    ! Propriedades físicas e Convectivas
    real(8), parameter :: Lx = 0.02, Ly = 0.02
    real(8), parameter :: rho = 2330.0, cp = 700.0, k = 150.0
    real(8), parameter :: alpha = k/(rho*cp)
    real(8), parameter :: h = 10.0, Tinf = 300.0
    
    ! Hotspot
    real(8), parameter :: Q0 = 1.0e8, x0 = Lx/2.0, y0 = Ly/2.0, sigma = 0.002, periodo = 0.5

    ! Parâmetros da malha e tempo
    integer, parameter :: Nx = 100, Ny = 100
    real(8) :: dx, dy, dt, tmax
    integer :: Nt
    
    ! Matrizes
    real(8) :: T(Nx, Ny), T_new(Nx, Ny), Q_gauss(Nx, Ny)

    ! Variáveis auxiliares
    integer :: i, j, n, iter, total_iter = 0
    real(8) :: rx, ry, x, y, tempo, Q, S

    ! Variáveis Gauss-Seidel
    integer, parameter :: max_iter = 5000
    real(8) :: erro_norma, soma_erro, T_old_iter, tol, den

    ! Variáveis para análises
    real(8) :: Tmax_atual, t_inicio_cpu, t_fim_cpu, tempo_cpu
    logical :: atingiu_limiar = .false.
    real(8) :: tempo_limiar
    real(8), parameter :: temp_alvo = 302.0

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

    tol = 1.0e-5 
    den = 1.0 + 2.0*rx + 2.0*ry                 ! Denominador da equação do Backward Euler

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
    write(*,'(A)') "        SIMULADOR TERMICO 2D - BACKWARD EULER          "
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

    ! Salva o histórico no tempo
    open(20, file = "serie_tempo_btcs.dat", status="replace")

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

        T_new = T                              ! Chute inicial

        ! --- SOLVER ITERATIVO DE GAUSS-SEIDEL ---
        do iter = 1, max_iter
            soma_erro = 0.0
            
            ! Cálculo do BTCS nos nós internos
            do j = 2, Ny-1
                do i = 2, Nx-1

                    Q = Q_gauss(i,j) * S
                    T_old_iter = T_new(i,j)

                    ! Equação do Backward Euler
                    T_new(i, j) = (T(i, j) &
                    + rx * (T_new(i+1, j) + T_new(i-1, j)) &
                    + ry * (T_new(i, j+1) + T_new(i, j-1)) &
                    + (dt * Q) / (rho*cp)) / den

                    ! Norma relativa (Acumula o quadrado do erro)
                    soma_erro = soma_erro + ((T_new(i,j) - T_old_iter) / T_new(i,j))**2

                end do
            end do

            ! Calcula a raiz quadrada da soma total
            erro_norma = sqrt(soma_erro)

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
            
            ! Critério de parada do solver
            if (erro_norma < tol) exit

        end do

        if (iter >= max_iter) print *, "AVISO: Não convergiu no tempo ", tempo
        total_iter = total_iter + iter

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
    write(*,'(A, F10.3)')           " Media de Iteracoes/Passo : ", real(total_iter)/real(Nt)
    write(*,'(A, F10.2, A)')        " Temperatura Max Final    : ", maxval(T), " K"
    write(*,'(A, ES10.2)')          " Erro L2 (Norma) Final    : ", erro_norma
    write(*,'(A)') "======================================================="
    write(*,'(A)') " --> Sucesso: 'btcs.dat' e 'serie_tempo_btcs.dat' gerados."

    ! -------------------------------------------------------------------
    ! SALVAMENTO DOS RESULTADOS
    ! -------------------------------------------------------------------
    open(10, file = "backward_euler.dat", status="replace")

    do i = 1, Nx
        do j = 1, Ny
            write(10,*) (i-1)*dx, (j-1)*dy, T(i,j)
        end do
    end do

    close(10)
    close(20)

end program simulacao_termica_btcs