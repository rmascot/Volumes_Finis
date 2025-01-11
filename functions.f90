module functions

    use const
    use visualisation

    implicit none

    contains

    !################# Initialisation de U #################
    subroutine initialize(U, dx, g, phi_g)
        real(PR), dimension(:,:), intent(out) :: U
        real(PR), dimension(:), intent(out) :: phi_g
        real(PR), intent(in) :: dx, g
        integer :: i
        real(PR) :: x

        U = 1._PR
        do i = 2, size(U)/3 - 1
            x = (i - 1 - 0.5) * dx
            U(i, 1) = (1._PR - ((gamma - 1._PR)/(gamma)) * g*x)**(1._PR/(gamma - 1._PR))  ! Densité
            U(i, 2) = 0._PR  ! Quantité de mouvement
            U(i, 3) = (U(i, 1)**gamma) / (gamma - 1._PR)  ! Énergie totale
            phi_g(i) = g * x
        end do
        phi_g(1) = phi_g(2)
        phi_g(size(U)/3) = phi_g(size(U)/3-1) - dx*g
    end subroutine initialize

    !################# Erreur en norme infinie #################
    subroutine calculate_infinity_norm(U, U_interp, Us, nx, erreur_inf_P, erreur_inf_u)
        real(PR), dimension(:,:), intent(in) :: U, U_interp, Us
        integer, intent(in) :: nx
        real(PR), intent(out) :: erreur_inf_P, erreur_inf_u
        real(PR), dimension(:), allocatable :: erreur_P, erreur_u
        real(PR) :: delta_P, delta_u
        integer :: i

        allocate(erreur_P(nx), erreur_u(nx))

        erreur_inf_P = 0._PR
        erreur_inf_u = 0._PR

        do i = 2, nx+1
            delta_P = abs((gamma - 1._PR) * (U(i, 3) - 0.5 * U(i, 2)**2 / U(i, 1)) - Us(i, 1)**gamma - U_interp(i, 4))
            delta_u = abs(U(i, 2) / U(i, 1) - U_interp(i, 2)/U_interp(i, 1))

            if (delta_P > erreur_inf_P) then
                erreur_inf_P = delta_P
                print*, i-1, (gamma - 1._PR) * (U(i, 3) - 0.5 * U(i, 2)**2 / U(i, 1)) - Us(i, 1)**gamma, U_interp(i, 4)
            end if
            if (delta_u > erreur_inf_u) then
                erreur_inf_u = delta_u
                !print*, i-1, U(i, 2) / U(i, 1), U_interp(i, 2)/U_interp(i, 1)
            end if
        end do
        deallocate(erreur_P, erreur_u)
    end subroutine calculate_infinity_norm

    !################# Erreur en norme L2 #################
    subroutine calculate_L2_norm(U, U_interp, Us, nx, erreur_L2_P, erreur_L2_u, erreur_rho)
        real(PR), dimension(:,:), intent(in) :: U, U_interp, Us
        integer, intent(in) :: nx
        real(PR), intent(out) :: erreur_L2_P, erreur_L2_u, erreur_rho
        real(PR) :: sum_delta_P, sum_delta_u, sum_erreur_rho
        real(PR) :: delta_P, delta_u
        integer :: i
    
        sum_delta_P = 0._PR
        sum_delta_u = 0._PR
        sum_erreur_rho = 0._PR
    
        ! Calcul des erreurs
        do i = 2, nx+1
            delta_P = abs((gamma - 1._PR) * (U(i, 3) - 0.5 * U(i, 2)**2 / U(i, 1)) - Us(i, 1)**gamma - U_interp(i, 4))
            delta_u = abs(U(i, 2) / U(i, 1) - U_interp(i, 2) / U_interp(i, 1))
    
            sum_delta_P = sum_delta_P + delta_P**2
            sum_delta_u = sum_delta_u + delta_u**2
            sum_erreur_rho = sum_erreur_rho + (U(i, 1) - U_interp(i, 1))**2
        end do
    
        ! Calcul des normes L2
        erreur_L2_P = sqrt(sum_delta_P / nx)
        erreur_L2_u = sqrt(sum_delta_u / nx)
        erreur_rho = sqrt(sum_erreur_rho / nx)
    end subroutine calculate_L2_norm

    !################# Erreur en norme L1 #################
    subroutine calculate_L1_norm(U, U_interp, Us, nx, erreur_P, erreur_u, erreur_rho)
        real(PR), dimension(:,:), intent(in) :: U, U_interp, Us
        integer, intent(in) :: nx
        real(PR), intent(out) :: erreur_P, erreur_u, erreur_rho
        real(PR) :: sum_delta_P, sum_delta_u, sum_erreur_rho
        real(PR) :: delta_P, delta_u
        integer :: i
    
        sum_delta_P = 0._PR
        sum_delta_u = 0._PR
        sum_erreur_rho = 0._PR
    
        ! Calcul des erreurs
        do i = 2, nx+1
            delta_P = abs((gamma - 1._PR) * (U(i, 3) - 0.5 * U(i, 2)**2 / U(i, 1)) - Us(i, 1)**gamma - U_interp(i, 4))
            delta_u = abs(U(i, 2) / U(i, 1) - U_interp(i, 2) / U_interp(i, 1))

            sum_delta_P = sum_delta_P + delta_P
            !print*, i, (gamma - 1._PR) * (U(i, 3) - 0.5 * U(i, 2)**2 / U(i, 1)) - Us(i, 1)**gamma, U_interp(i, 4)
            sum_delta_u = sum_delta_u + delta_u
            sum_erreur_rho = sum_erreur_rho + abs(U(i, 1) - U_interp(i, 1))
        end do
    
        ! Calcul des normes L2
        erreur_P = sum_delta_P / nx
        erreur_u = sum_delta_u / nx
        erreur_rho = sum_erreur_rho / nx
    end subroutine calculate_L1_norm

    !################# Initialisation des cellules fantômes pour les conditions de bord #################
    subroutine set_ghost_cells(U, P, nx)
        real(PR), dimension(:,:), intent(inout) :: U
        real(PR), dimension(:), intent(inout) :: P
        integer, intent(in) :: nx
    
        ! Cellule fantôme à gauche
        U(1, :) = U(2, :)  ! Extrapolation de densité, quantité de mouvement, énergie
        P(1) = P(2)        ! Extrapolation de la pression
    
        ! Cellule fantôme à droite
        U(nx+2, :) = U(nx+1, :)
        P(nx+2) = P(nx+1)

    end subroutine set_ghost_cells

    !################# Calcul de - div (F) #################
    subroutine div_flux(F_arete, m_div_F, dx)
        !real(PR), dimension(:), intent(in) :: Flux_gauche, Flux_droite
        real(PR), dimension(:,:), intent(in) :: F_arete
        real(PR), dimension(:,:), intent(inout) :: m_div_F
        real(PR), intent(in) :: dx
        integer :: i
    
        ! Calcul des flux de Rusanov
        do i = 1, size(F_arete)/3 - 1
            ! Ajout du - div(F)
            m_div_F(i, :) = - ( F_arete(i+1, :) - F_arete(i, :) ) / (dx)
        end do

    end subroutine div_flux

    !################# Calcul du flux sur une maille #################
    subroutine calcul_flux_arete(U, F, F_arete, max_vp)
        real(PR), dimension(:,:), intent(in) :: U
        real(PR), dimension(:), intent(in) :: max_vp
        real(PR), dimension(:,:), intent(in) :: F
        real(PR), dimension(:,:), intent(inout) :: F_arete
        integer :: i

        do i = 1, size(F_arete)/3
            F_arete(i, :) = 0.5*(F(i+1, :) + F(i, :)) - 0.5 * max(max_vp(i), max_vp(i+1)) * (U(i+1, :) - U(i, :))
        end do

    end subroutine calcul_flux_arete

    !################# Calcul du flux sur les arêtes #################
    subroutine calcul_flux_maille (U, P, F)
        real(PR), dimension(:,:), intent(in) :: U
        real(PR), dimension(:), intent(in) :: P
        real(PR), dimension(:,:), intent(inout) :: F
        integer :: i

        do i = 1, size(F)/3
            F(i, 1) = U(i, 2)                               ! Flux de masse
            F(i, 2) = (U(i, 2)**2 / U(i, 1)) + P(i)         ! Flux de quantité de mouvement
            F(i, 3) = U(i, 2)/U(i, 1) * ((U(i, 3) + P(i)))  ! Flux d'énergie
        end do

    end subroutine calcul_flux_maille

    !################# Calcul du Terme Source #################
    subroutine compute_source(U, source, g)
        real(PR), dimension(:,:), intent(in) :: U
        real(PR), dimension(:,:), intent(out) :: source
        real(PR), intent(in) :: g
        integer :: i
    
        do i = 1, size(source)/3
            source(i, 1) = 0.0     
            source(i, 2) = - U(i+1, 1) * g   
            source(i, 3) = - U(i+1, 2) * g 
        end do
    end subroutine compute_source

    !################# Lecture des parametres #################
    subroutine read_parameters(filename, nx, g, L, tmax)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: nx
        real(PR), intent(out) :: g, L, tmax
        character(len=100) :: line
        integer :: ios
    
        open(unit=10, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            print*, "Error opening file: ", filename
            stop
        end if
    
        do
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) exit
    
            select case (adjustl(trim(line)))
            case ("nx =")
                read(10, *) nx
            case ("g =")
                read(10, *) g
            case ("L =")
                read(10, *) L
            case ("tmax =")
                read(10, *) tmax
            end select
        end do
    
        close(10)
    end subroutine

end module functions