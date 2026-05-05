! =============================================================
! Simulador térmico 2D - LAPLACE
! =============================================================

program simulador_laplace
    implicit none

    ! Parâmetros do domínio
    real(8) :: Lx, Ly, tol, dx, dy, aE, aW, aN, aS, aP
    integer :: Nx, Ny, max_iter
    character(len=*), parameter :: out_dir = '../data/results/'

    ! Variáveis do solver
    real(8) :: T_new, R, residuo_max
    integer :: i, j, iter

    real(8), allocatable :: T(:,:)

    call init_parametros()
    allocate(T(Nx, Ny))

    T = 0.0d0

    ! Condições de Contorno de Dirichlet
    T(:, Ny) = 100.0d0 ! y = Ly (Topo)
    T(:, 1)  = 0.0d0   ! y = 0 (Base)
    T(1, :)  = 0.0d0   ! x = 0 (Esquerda)
    T(Nx, :) = 0.0d0   ! x = Lx (Direita)

    open(10, file=trim(out_dir)//'norma_iter_laplace.dat', status='replace')

    ! Gauss-Seidel
    do iter = 1, max_iter
        residuo_max = 0.0d0

        do j = 2, Ny-1
            do i = 2, Nx-1
                T_new = (aE*T(i+1,j) + aW*T(i-1,j) + aN*T(i,j+1) + aS*T(i,j-1)) / aP

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

    ! Mapa de Calor
    open(50, file=trim(out_dir)//'mapa_calor_laplace.dat', status='replace')
    do j = 1, Ny
        do i = 1, Nx
            write(50, '(2F12.6, ES15.6)') (i-1)*dx, (j-1)*dy, T(i,j)
        end do
        write(50, *)
    end do
    close(50)

    print *, "Simulacao de Laplace Concluida!"
    deallocate(T)

contains

    subroutine init_parametros()
        namelist /parametros/ Lx, Ly, Nx, Ny, tol, max_iter

        open(99, file='../parametros_laplace.txt', status='old')
        read(99, nml=parametros)
        close(99)
        
        dx = Lx / dble(Nx - 1)
        dy = Ly / dble(Ny - 1)

        aE = dy**2
        aW = dy**2
        aN = dx**2
        aS = dx**2
        aP = 2.0d0 * (dx**2 + dy**2)

        call execute_command_line('mkdir -p ' // trim(out_dir))
    end subroutine init_parametros

end program simulador_laplace