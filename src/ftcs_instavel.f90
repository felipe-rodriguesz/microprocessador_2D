program ftcs_instavel
    implicit none

    ! -------------------------------------------------------------------
    ! DECLARAÇÕES DE VARIÁVEIS E CONSTANTES
    ! -------------------------------------------------------------------
    real(8), parameter :: Lx = 0.02, Ly = 0.02
    real(8), parameter :: rho = 2330.0, cp = 700.0, k = 150.0
    real(8), parameter :: alpha = k/(rho*cp)
    real(8), parameter :: h = 10.0, Tinf = 300.0
    real(8), parameter :: Q0 = 1.0e8, x0 = Lx/2.0, y0 = Ly/2.0
    real(8), parameter :: sigma = 0.002, periodo = 0.5

    integer, parameter :: Nx = 100, Ny = 100
    real(8) :: dx, dy, dt, dt_estavel, tmax
    integer :: Nt

    real(8) :: T(Nx, Ny), T_new(Nx, Ny), Q_gauss(Nx, Ny)
    integer :: i, j, n
    real(8) :: rx, ry, x, y, tempo, Q, S
    real(8) :: Tmax_atual
    logical :: divergiu = .false.
    real(8), parameter :: T_limite = 1.0e6   ! Limiar de detecção de divergência (K)

    ! -------------------------------------------------------------------
    ! CONFIGURAÇÃO DA MALHA
    ! -------------------------------------------------------------------
    dx = Lx / (Nx - 1)
    dy = Ly / (Ny - 1)

    ! dt estável original
    dt_estavel = (dx**2) / (4.0 * alpha)

    ! dt INSTÁVEL — 3x maior que o limite
    dt = dt_estavel * 3.0

    tmax = 4.75
    Nt = int(tmax / dt)

    rx = alpha * dt / (dx**2)
    ry = alpha * dt / (dy**2)

    ! -------------------------------------------------------------------
    ! IMPRESSÃO DO AVISO
    ! -------------------------------------------------------------------
    write(*,'(A)') "========================================================="
    write(*,'(A)') "   FTCS INSTÁVEL — APENAS PARA DEMONSTRAÇÃO DIDÁTICA"
    write(*,'(A)') "========================================================="
    write(*,'(A, F8.6, A)') " dt estável original : ", dt_estavel, " s"
    write(*,'(A, F8.6, A)') " dt instável usado   : ", dt,         " s  (3x maior)"
    write(*,'(A, F8.4)')    " rx                  : ", rx
    write(*,'(A, F8.4)')    " ry                  : ", ry
    write(*,'(A, F8.4)')    " rx + ry             : ", rx + ry
    write(*,'(A)') ""
    write(*,'(A)') " Condicao de Von Neumann: rx + ry <= 0.50"
    if (rx + ry > 0.5) then
        write(*,'(A, F5.2, A)') " STATUS: VIOLADA (rx + ry = ", rx+ry, ") — instabilidade esperada!"
    end if
    write(*,'(A)') "========================================================="
    write(*,'(A)') " Iniciando simulação... aguarde divergência."
    write(*,'(A)') ""

    ! -------------------------------------------------------------------
    ! CONDIÇÃO INICIAL E PRÉ-CÁLCULO
    ! -------------------------------------------------------------------
    T = Tinf
    T_new = Tinf
    T(:,1) = Tinf

    do j = 1, Ny
        y = (j-1)*dy
        do i = 1, Nx
            x = (i-1)*dx
            Q_gauss(i,j) = Q0 * exp(-((x-x0)**2 + (y-y0)**2)/(2.0*sigma**2))
        end do
    end do

    ! -------------------------------------------------------------------
    ! SALVAMENTO DA SÉRIE TEMPORAL
    ! -------------------------------------------------------------------
    open(20, file = "serie_tempo_ftcs_instavel.dat", status="replace")

    ! ===================================================================
    ! LOOP TEMPORAL
    ! ===================================================================
    do n = 1, Nt
        tempo = dble(n) * dt

        ! Duty-Cycle por tempo físico
        if ( mod(tempo, periodo) < (periodo / 2.0d0) ) then
            S = 1.0
        else
            S = 0.0
        end if

        ! FTCS nos nós internos
        do j = 2, Ny-1
            do i = 2, Nx-1
                Q = Q_gauss(i,j) * S

                T_new(i,j) = T(i,j) &
                + rx * (T(i+1,j) - 2.0*T(i,j) + T(i-1,j)) &
                + ry * (T(i,j+1) - 2.0*T(i,j) + T(i,j-1)) &
                + dt * Q / (rho*cp)
            end do
        end do

        ! Condições de contorno
        do j = 1, Ny
            T_new(1,j)  = T_new(2,j)
            T_new(Nx,j) = T_new(Nx-1,j)
        end do
        T_new(:,1) = Tinf
        do i = 1, Nx
            T_new(i,Ny) = (k/dy*T_new(i,Ny-1) + h*Tinf)/(k/dy + h)
        end do

        T = T_new

        Tmax_atual = maxval(abs(T))

        ! Salva a série temporal
        write(20, *) tempo, T(Nx/2, Ny/2), Tmax_atual

        ! Detecta divergência e para antes de gerar NaN/Inf
        if (Tmax_atual > T_limite .or. isnan(Tmax_atual)) then
            write(*,'(A, I6, A)')      " DIVERGÊNCIA detectada no passo n = ", n, "!"
            write(*,'(A, F10.4, A)')   " Tempo                            = ", tempo, " s"
            write(*,'(A, ES12.4, A)')  " T_max atingiu                    = ", Tmax_atual, " K"
            write(*,'(A)') " Simulação interrompida para evitar overflow."
            divergiu = .true.
            exit
        end if

    end do

    close(20)

    ! -------------------------------------------------------------------
    ! RELATÓRIO FINAL
    ! -------------------------------------------------------------------
    write(*,'(A)') ""
    write(*,'(A)') "========================================================="
    write(*,'(A)') "              RELATÓRIO — FTCS INSTÁVEL"
    write(*,'(A)') "========================================================="
    if (divergiu) then
        write(*,'(A)') " RESULTADO: Instabilidade numérica confirmada."
        write(*,'(A)') " A temperatura divergiu exponencialmente,"
        write(*,'(A)') " comprovando a violação da condição de Von Neumann."
    else
        write(*,'(A)') " RESULTADO: Não divergiu no tempo simulado."
        write(*,'(A)') " Tente um dt ainda maior para forçar a instabilidade."
    end if
    write(*,'(A)') "========================================================="
    write(*,'(A)') " --> Arquivo 'serie_tempo_ftcs_instavel.dat' gerado."

    ! Campo final (pode conter valores gigantes — salvo para visualização)
    open(10, file = "ftcs_instavel.dat", status="replace")
    do i = 1, Nx
        do j = 1, Ny
            write(10,*) (i-1)*dx, (j-1)*dy, T(i,j)
        end do
    end do
    close(10)

end program ftcs_instavel
