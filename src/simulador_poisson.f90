! =============================================================
! Simulador térmico 2D - POISSON
! Objetivo: MMS com T(x,y)=sin(pi x) cos(pi y)
! =============================================================

program simulador_poisson
    implicit none

    ! Parâmetros de entrada e derivados
    real(8) :: Lx, Ly, tol, dx, dy, aE, aW, aN, aS, aP
    integer :: Nx, Ny, max_iter
    character(len=*), parameter :: out_dir = '../data/results/'

    ! Variáveis do solver
    real(8) :: T_new, R, residuo_max
    real(8) :: pos_x, pos_y, T_exato, erro_atual, erro_max
    real(8) :: kx, ky
    real(8), parameter :: pi = acos(-1.0d0)
    integer :: i, j, iter

    real(8), allocatable :: T(:,:), b(:,:)

    call init_parametros()
    allocate(T(Nx, Ny), b(Nx, Ny))

    T = 0.0d0  ! Chute inicial uniforme

    ! Arquivo de convergencia
    open(10, file=trim(out_dir)//'norma_iter_poisson.dat', status='replace')

    ! Parametros da solucao exata
    kx = pi / Lx
    ky = pi / Ly

    ! Termo fonte 
    do j = 2, Ny-1
        do i = 2, Nx-1
            pos_x = (i-1)*dx
            pos_y = (j-1)*dy
            b(i,j) = -(kx**2 + ky**2) * sin(kx*pos_x) * cos(ky*pos_y)
        end do
    end do

    ! Condições de Contorno (Dirichlet)
    do i = 1, Nx
        pos_x = (i-1)*dx
        T(i,1)  = sin(kx*pos_x)
        T(i,Ny) = sin(kx*pos_x) * cos(ky*Ly)
    end do
    do j = 1, Ny
        pos_y = (j-1)*dy
        T(1,j)  = 0.0d0
        T(Nx,j) = sin(kx*Lx) * cos(ky*pos_y)
    end do

    ! Gauss-Seidel
    do iter = 1, max_iter
        residuo_max = 0.0d0

        do j = 2, Ny-1
            do i = 2, Nx-1
                T_new = (aE*T(i+1,j) + aW*T(i-1,j) + aN*T(i,j+1) + aS*T(i,j-1) - b(i,j)) / aP
                R = abs(T_new - T(i,j))
                if (R > residuo_max) residuo_max = R
                T(i,j) = T_new
            end do
        end do

        write(10,*) iter, residuo_max
        if (residuo_max < tol) exit
    end do

    write(*,'(A,I5,A,ES10.3)') "Iteracoes GS: ", iter, " | Residuo Maximo: ", residuo_max

    close(10)

    ! Erro maximo contra a solucao exata
    erro_max = 0.0d0
    do j = 1, Ny
        do i = 1, Nx
            pos_x = (i-1)*dx
            pos_y = (j-1)*dy
            T_exato = sin(kx*pos_x) * cos(ky*pos_y)
            erro_atual = abs(T(i,j) - T_exato)
            if (erro_atual > erro_max) erro_max = erro_atual
        end do
    end do

    ! Mapa de Calor
    open(50, file=trim(out_dir)//'mapa_calor_poisson.dat', status='replace')
    do j = 1, Ny
        do i = 1, Nx
            write(50, '(2F12.6, ES15.6)') (i-1)*dx, (j-1)*dy, T(i,j)
        end do
    end do
    close(50)

    print *, "==================================================="
    print *, "Teste de Poisson Concluido!"
    print *, "Erro Maximo contra a solucao exata: ", erro_max
    print *, "==================================================="
    
    deallocate(T, b)

contains

    subroutine init_parametros()
        namelist /parametros/ Lx, Ly, Nx, Ny, tol, max_iter

        open(99, file='../parametros_poisson.txt', status='old')
        read(99, nml=parametros)
        close(99)
        
        ! Calculos derivados e coeficientes
        dx = Lx / dble(Nx - 1)
        dy = Ly / dble(Ny - 1)

        aE = 1.0d0 / dx**2
        aW = 1.0d0 / dx**2
        aN = 1.0d0 / dy**2
        aS = 1.0d0 / dy**2
        aP = aE + aW + aN + aS

        call execute_command_line('mkdir -p ' // trim(out_dir))
    end subroutine init_parametros

end program simulador_poisson