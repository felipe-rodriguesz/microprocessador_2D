program crank_nicolson
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
    real(8) :: erro_norma, soma_erro, T_old_iter, tol, den
    real(8) :: residuo_max, residuo_local, vizinhos

    ! Variáveis para análises
    real(8) :: Tmax_atual, t_inicio_cpu, t_fim_cpu, tempo_cpu
    logical :: atingiu_limiar = .false.
    real(8) :: tempo_limiar
    real(8), parameter :: temp_alvo = 302.0

    ! Variáveis para os Snapshots
    integer :: proximo_segundo = 1
    character(len=30) :: nome_arquivo

    ! -------------------------------------------------------------------
    ! ÁREA DE EXECUÇÃO (Cálculos e Loops)
    ! -------------------------------------------------------------------

    ! Pede o tamanho da malha
    write(*,*) "======================================================="
    write(*,*) " Digite o tamanho da malha espacial (Nx e Ny):"
    write(*,*) " Exemplo para 100x100, digite: 100 100"
    write(*,*) "======================================================="
    read(*,*) Nx, Ny, fator_dt
    allocate( T(Nx, Ny), T_new(Nx, Ny), Q_gauss(Nx, Ny), RHS(Nx, Ny) )

    dx = Lx / dble(Nx - 1)
    dy = Ly / dble(Ny - 1)
    
    dt = fator_dt * ((dx**2) / (4.0 * alpha))
    tmax = 4.75
    Nt = int(tmax / dt)

    rx = alpha * dt / (dx**2) 
    ry = alpha * dt / (dy**2)

    ! Fatores do Crank-Nicolson (metade do FTCS)
    Fx = rx / 2.0
    Fy = ry / 2.0

    tol = 1.0e-5 
    den = 1.0 + 2.0*Fx + 2.0*Fy                 ! Denominador da equação do Crank-Nicolson

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
    write(*,'(A)') "======================================================="
    write(*,'(A)') " Iniciando loop temporal... aguarde."
    write(*,'(A)') ""

    ! Salva o histórico no tempo
    open(20, file = "serie_tempo_cn.dat", status="replace")

    ! Cronômetro do processador
    call cpu_time(t_inicio_cpu)

    ! ===================================================================
    ! Loop temporal
    ! ===================================================================
    do n = 1, Nt
        tempo = n * dt

        ! Estado atual do Duty-Cycle (instante t)
        if ( mod(tempo, periodo) < (periodo / 2.0) ) then
            S = 1.0
        else
            S = 0.0
        end if

        ! Estado futuro do Duty-Cycle (instante t + dt)
        if ( mod(tempo + dt, periodo) < (periodo / 2.0) ) then
            S_next = 1.0
        else
            S_next = 0.0
        end if

        ! CALCULAR O LADO DIREITO (RHS - Parte Explícita)
        ! Usa as temperaturas atuais (T) para calcular a inércia térmica
        do j = 2, Ny-1
            do i = 2, Nx-1
                Q = Q_gauss(i,j) * S
                Q_next = Q_gauss(i,j) * S_next
                
                RHS(i,j) = T(i,j) &
                         + Fx * (T(i+1,j) - 2.0*T(i,j) + T(i-1,j)) &
                         + Fy * (T(i,j+1) - 2.0*T(i,j) + T(i,j-1)) &
                         + (dt * (Q + Q_next)) / (2.0 * rho * cp)
            end do
        end do

        T_new = T                              ! Chute inicial

        ! --- SOLVER ITERATIVO DE GAUSS-SEIDEL 
        do iter = 1, max_iter
            soma_erro = 0.0
            residuo_max = 0.0
            
            ! Cálculo do Crank-Nicolson nos nós internos
            do j = 2, Ny-1
                do i = 2, Nx-1

                    T_old_iter = T_new(i,j)

                    ! Soma as influências dos vizinhos recém-calculados
                    vizinhos = Fx * (T_new(i+1, j) + T_new(i-1, j)) &
                             + Fy * (T_new(i, j+1) + T_new(i, j-1))

                    ! Calcula o Resíduo de Gauss-Seidel (R = b - Ax)
                    residuo_local = RHS(i,j) + vizinhos - (den * T_old_iter)

                    ! Atualiza a temperatura local usando o resíduo
                    T_new(i, j) = T_old_iter + (residuo_local / den)

                    ! Acumula os erros para análise global
                    soma_erro = soma_erro + residuo_local**2
                    if (abs(residuo_local) > residuo_max) residuo_max = abs(residuo_local)

                end do
            end do

            ! Calcula o Resíduo RMS
            erro_norma = sqrt(soma_erro / dble(Nx * Ny))

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

        ! Salvar campos intermediários (Snapshots)
        ! Verifica se o tempo cruzou 1.0, 2.0, 3.0 ou 4.0 segundos
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
    write(*,'(A, ES10.2)')          " Residuo Maximo Final     : ", residuo_max
    write(*,'(A, ES10.2)')          " Erro L2 (Norma) Final    : ", erro_norma
    write(*,'(A)') "======================================================="
    write(*,'(A)') " --> Sucesso: 'crank_nicolson.dat' e 'serie_tempo_cn.dat' gerados."

    ! -------------------------------------------------------------------
    ! SALVAMENTO DOS RESULTADOS
    ! -------------------------------------------------------------------
    open(10, file = "crank_nicolson.dat", status="replace")

    do i = 1, Nx
        do j = 1, Ny
            write(10,*) (i-1)*dx, (j-1)*dy, T(i,j)
        end do
    end do

    close(10)
    close(20)

    ! -------------------------------------------------------------------
    ! EXPORTAÇÃO PARA ESTUDO DE CONVERGÊNCIA
    ! -------------------------------------------------------------------
    open(40, file = "dados_convergencia.dat", position="append", status="unknown")
    
    ! Salva: [1] Número de Nós, [2] Tamanho do dx, [3] Temperatura Máxima
    write(40, *) Nx, dx, maxval(T)
    
    close(40)
    print *, "--> Ponto de convergencia salvo em 'dados_convergencia.dat'."

    ! -------------------------------------------------------------------
    ! SALVAMENTO DE MÉTRICAS E LOGS
    ! -------------------------------------------------------------------
    ! LOG ESTRUTURADO (.txt)
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
    write(50,'(A, F10.2, A)')     " T_max final      : ", maxval(T), " K"
    write(50,'(A, F10.3)')        " Media Iter/Passo : ", real(total_iter)/real(Nt)
    write(50,'(A, ES10.2)')       " Residuo Maximo   : ", residuo_max
    write(50,'(A, ES10.2)')       " Erro L2 (Norma)  : ", erro_norma
    if (atingiu_limiar) then
        write(50,'(A, F10.3, A)') " Limiar (302K)    : ", tempo_limiar, " s"
    else
        write(50,'(A)')           " Limiar (302K)    : Nao atingido"
    end if
    write(50,'(A)')         "======================================================="
    close(50)

    ! DADOS BRUTOS PARA O MATLAB (.dat)
    open(51, file = "metricas_cn.dat", status="replace")
    if (atingiu_limiar) then
        write(51, *) Nx, dx, dt, rx, tempo_cpu, maxval(T), tempo_limiar, real(total_iter)/real(Nt), residuo_max, erro_norma
    else
        write(51, *) Nx, dx, dt, rx, tempo_cpu, maxval(T), -1.0d0, real(total_iter)/real(Nt), residuo_max, erro_norma
    end if
    close(51)

    ! Libera a memória RAM utilizada pelas matrizes
    deallocate( T, T_new, Q_gauss, RHS )

end program crank_nicolson