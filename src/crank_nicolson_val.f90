program crank_nicolson_val
    implicit none

    ! Variáveis Estritamente Necessárias para a Matemática
    real(8) :: Lx, Ly, rho, cp, k, alpha, tmax
    integer :: Nx, Ny, Nt
    real(8) :: dx, dy, dt, fator_dt

    ! Matrizes
    real(8), allocatable :: T(:, :), T_new(:, :), RHS(:, :)

    ! Variáveis auxiliares
    integer :: i, j, n, iter
    real(8) :: rx, ry, x, y, tempo, Fx, Fy, pi
    real(8) :: T_exato, erro_max

    ! Variáveis Gauss-Seidel
    integer, parameter :: max_iter = 5000
    real(8) :: T_old_iter, tol, den, residuo_local, vizinhos, norma

    ! ===================================================================
    ! LEITURA DO ARQUIVO DE TESTE
    ! ===================================================================
    open(15, file="parametros_teste.txt", status="old", action="read")
    read(15, *) Lx, Ly
    read(15, *) rho, cp, k
    read(15, *) tmax
    read(15, *) Nx, Ny, fator_dt
    close(15)

    ! Pi e Difusividade
    pi = 4.0d0 * atan(1.0d0)
    alpha = k / (rho * cp)
    
    ! Discretização Espacial
    dx = Lx / dble(Nx - 1)
    dy = Ly / dble(Ny - 1)
    
    ! Discretização Temporal
    dt = fator_dt * ((dx**2) / (4.0d0 * alpha))
    Nt = int(tmax / dt) + 2  ! +2 para garantir que o tempo cruza o tmax

    rx = alpha * dt / (dx**2) 
    ry = alpha * dt / (dy**2)

    Fx  = rx / 2.0d0
    Fy  = ry / 2.0d0
    tol = 1.0e-6 
    den = 1.0d0 + 2.0d0*Fx + 2.0d0*Fy

    allocate( T(Nx, Ny), T_new(Nx, Ny), RHS(Nx, Ny) )

    ! ===================================================================
    ! CONDIÇÃO INICIAL (Lombada Senoidal)
    ! ===================================================================
    do j = 1, Ny
        y = (j-1)*dy
        do i = 1, Nx
            x = (i-1)*dx
            T(i,j) = sin(pi * x / Lx) * sin(pi * y / Ly)
        end do
    end do
    T_new = T

    ! ===================================================================
    ! LOOP TEMPORAL (Crank-Nicolson)
    ! ===================================================================
    do n = 1, Nt
        tempo = n * dt

        ! Lado direito (RHS) sem fonte de calor
        do j = 2, Ny-1
            do i = 2, Nx-1
                RHS(i,j) = T(i,j) &
                         + Fx * (T(i+1,j) - 2.0d0*T(i,j) + T(i-1,j)) &
                         + Fy * (T(i,j+1) - 2.0d0*T(i,j) + T(i,j-1))
            end do
        end do

        ! Solver Gauss-Seidel
        do iter = 1, max_iter
            norma = 0.0d0

            do j = 2, Ny-1
                do i = 2, Nx-1
                    T_old_iter = T_new(i,j)
                    vizinhos = Fx * (T_new(i+1, j) + T_new(i-1, j)) &
                             + Fy * (T_new(i, j+1) + T_new(i, j-1))

                    residuo_local = RHS(i,j) + vizinhos - (den * T_old_iter)
                    T_new(i, j) = T_old_iter + (residuo_local / den)
                    norma = max(norma, abs(T_new(i,j) - T_old_iter))
                end do
            end do

            ! Contornos de Dirichlet (T = 0 nas bordas)
            T_new(1,:) = 0.0d0
            T_new(Nx,:) = 0.0d0
            T_new(:,1) = 0.0d0
            T_new(:,Ny) = 0.0d0

            if (norma < tol) exit
        end do

        T = T_new
        if (tempo >= tmax) exit
    end do

    ! ===================================================================
    ! CÁLCULO DO ERRO MÁXIMO E GRAVAÇÃO
    ! ===================================================================
    erro_max = 0.0d0
    do j = 1, Ny
        y = (j-1)*dy
        do i = 1, Nx
            x = (i-1)*dx
            ! Equação Analítica Exata
            T_exato = sin(pi * x / Lx) * sin(pi * y / Ly) * &
                      exp(-alpha * (pi**2) * (1.0d0/(Lx**2) + 1.0d0/(Ly**2)) * tempo)
            
            erro_max = max(erro_max, abs(T(i,j) - T_exato))
        end do
    end do

    ! Salva os resultados
    open(70, file = "erro_convergencia.dat", position="append", status="unknown")
    write(70, '(I5, 3ES15.6)') Nx, dx, dt, erro_max
    close(70)

    print *, "========================================="
    print *, " VALIDAÇÃO ANALÍTICA - CRANK-NICOLSON"
    print *, "========================================="
    print *, "Malha (Nx) : ", Nx
    print *, "Delta x    : ", dx
    print *, "Delta t    : ", dt
    print *, "Erro Maximo: ", erro_max
    print *, "========================================="

    deallocate( T, T_new, RHS )

end program crank_nicolson_val