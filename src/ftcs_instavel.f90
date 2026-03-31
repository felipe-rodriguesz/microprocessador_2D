program ftcs_instavel
    implicit none

    ! -------------------------------------------------------------------
    ! DECLARAÇÕES DE VARIÁVEIS E CONSTANTES
    ! -------------------------------------------------------------------
    ! Variáveis Físicas
    real(8) :: Lx, Ly
    real(8) :: rho, cp, k, alpha
    real(8) :: h, Tinf
    real(8) :: sigma, periodo, tmax

    ! Hotspot e malha
    real(8) :: Q0, x0, y0, duty, fator_dt
    integer :: Nx, Ny
    real(8) :: dx, dy, dt, dt_estavel
    integer :: Nt

    real(8), allocatable :: T(:, :), T_new(:, :), Q_gauss(:, :)
    integer :: i, j, n
    real(8) :: rx, ry, x, y, tempo, Q, S
    real(8) :: Tmax_atual
    logical :: divergiu = .false.
    real(8), parameter :: T_limite = 1.0e6   ! Limiar de detecção de divergência

    read(*,*) Nx, Ny, fator_dt, Q0, x0, y0, duty

    allocate(T(Nx, Ny), T_new(Nx, Ny), Q_gauss(Nx, Ny))

    ! -------------------------------------------------------------------
    ! CONFIGURAÇÃO DA MALHA E LEITURA DO ARQUIVO DE PARÂMETROS
    ! -------------------------------------------------------------------
    open(15, file="src/parametros.txt", status="old", action="read")
    read(15, *) Lx, Ly
    read(15, *) rho, cp, k
    read(15, *) h, Tinf
    read(15, *) sigma, periodo, tmax
    close(15)

    alpha = k / (rho * cp)

    dx = Lx / (Nx - 1)
    dy = Ly / (Ny - 1)

    ! dt estável original
    dt_estavel = (dx**2) / (4.0 * alpha)

    ! dt INSTÁVEL — 3x maior que o limite (preserva fator_dt via stdin)
    dt = dt_estavel * 3.0d0 * fator_dt

    Nt = int(tmax / dt)

    rx = alpha * dt / (dx**2)
    ry = alpha * dt / (dy**2)

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

        ! Duty-Cycle no instante t
        if ( mod(tempo, periodo) < ((1.0 - duty) * periodo) ) then
            S = 0.0  ! Começa DESLIGADA
        else
            S = 1.0  ! Termina LIGADA
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
            divergiu = .true.
            exit
        end if

    end do

    close(20)

    ! Campo final (pode conter valores gigantes — salvo para visualização)
    open(10, file = "ftcs_instavel.dat", status="replace")
    do i = 1, Nx
        do j = 1, Ny
            write(10,*) (i-1)*dx, (j-1)*dy, T(i,j)
        end do
    end do
    close(10)

    deallocate(T, T_new, Q_gauss)

end program ftcs_instavel