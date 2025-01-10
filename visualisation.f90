module visualisation

    use const

    implicit none

    contains

    !################# Écriture des résultats #################
    subroutine output_results(U, Us, P, nx, L)
        real(PR), dimension(:,:), intent(in) :: U, Us
        real(PR), dimension(:), intent(in) :: P
        real(PR), intent(in) :: L
        integer, intent(in) :: nx
        integer :: i
        real(PR) :: x
        character(len=20) :: filename
    
        filename = "results.txt"
        open(unit=10, file=filename, status="replace")
    
        write(10,*) "x, rho, rho*u, E, delta_P, delta_u, u_init, P_init"
        do i = 2, nx+1
            x = (i - 1 - 0.5) * L / nx
            write(10,*) x, U(i, 1), U(i, 2), U(i, 3), P(i)-Us(i,1)**gamma, (U(i, 2)/U(i, 1)), (Us(i, 2)/Us(i, 1)), Us(i,1)**gamma
            !write(10,*) x, Us(i, 1), Us(i, 2)/Us(i, 1), U(i, 3), Us(i,1)**gamma, (Us(i, 2)/Us(i, 1))
        end do
    
        close(10)
        print *, "Résultats enregistrés dans : ", filename
    end subroutine output_results

    !################# Calcul de la solution exacte #################
    subroutine compute_analytical(U_analytical, t, nx, dx, L)
        real(PR), dimension(:,:), intent(out) :: U_analytical
        real(PR), intent(in) :: t, dx, L
        integer, intent(in) :: nx
        integer :: i
        real(PR) :: x, rho_0, u_0, p_0, epsilon, c_s, k
    
        ! Paramètres de l'onde
        rho_0 = 1.0
        u_0 = 0.0
        p_0 = 1.0
        epsilon = 1.0
        c_s = sqrt(gamma * p_0 / rho_0)
        k = 2.0 * pi / L
    
        ! Calcul de la solution analytique
        do i = 1, nx
            x = (i - 0.5) * dx
            U_analytical(i, 1) = rho_0 + epsilon * sin(k * x - c_s * t)  ! Densité
            U_analytical(i, 2) = (rho_0 + epsilon * sin(k * x - c_s * t)) * &
                                 (u_0 + epsilon * c_s / rho_0 * sin(k * x - c_s * t))  ! Quantité de mouvement
            U_analytical(i, 3) = p_0 / (gamma - 1) + 0.5 * U_analytical(i, 2)**2 / U_analytical(i, 1)  ! Énergie totale
        end do
    end subroutine compute_analytical

end module visualisation