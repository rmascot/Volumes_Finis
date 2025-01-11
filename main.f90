program euler_rusanov_gravity

    ! Utilisation des modules
    use const
    use visualisation
    use functions

    implicit none

    ! Variables pour les calculs
    logical :: test_convergence
    real(PR) :: dx, L, g, tmax, erreur_P, erreur_u
    real(PR), dimension(:,:), allocatable :: U, U_ref, U_interp, Us
    real(PR), dimension(:), allocatable :: phi_g, P
    integer :: nx, nx_ref, i
    real(PR) :: temp1, temp2, temp3, temp4, temp5, temp6  ! Variables pour lire les colonnes

    ! Initialisation des tableaux pour stocker les erreurs
    real(PR), dimension(5) :: erreurs_P, erreurs_u
    integer, dimension(5) :: nx_values
    character(len=256) :: filename

    ! === PARAMÈTRES UTILISATEUR ===
    test_convergence = .TRUE.  ! Mettre à .FALSE. pour une seule simulation

    if (.not. test_convergence) then
        ! =============================
        ! CAS SIMPLE : UNE SEULE SIMULATION
        ! =============================

        call read_parameters("parameters.txt", nx, g, L, tmax)

        dx = L / nx
        allocate(U(nx+2, 3), phi_g(nx+2), P(nx+2))

        call initialize(U, dx, g, phi_g)
        call simulate(U, P, nx, g, L, tmax, phi_g)

        allocate(Us(nx+2, 3))
        call initialize(Us, dx, g, phi_g)
        call output_results(U, Us, P, nx, L, "solution_reference.txt")

        deallocate(U, P, phi_g, Us)
        print*, "Simulation terminée avec succès."

    else
        ! =============================
        ! CAS TEST DE CONVERGENCE
        ! =============================

        call read_parameters("parameters.txt", nx_ref, g, L, tmax)

        allocate(U_ref(nx_ref, 4))

        ! Ouvrir le fichier
        open(unit=10, file="solution_reference.txt", status="old", action="read")
    
        ! Lire les données
        do i = 1, nx_ref
            read(10, *) temp1, temp2, temp3, temp4, temp5, temp6
            U_ref(i, 1) = temp2                   ! rho
            U_ref(i, 2) = temp3                   ! rho*u
            U_ref(i, 3) = temp4                   ! E
            U_ref(i, 4) = temp5                   ! P
        end do

        do i = 1, 5
            nx = 128*2**(i-1)
            nx_values(i) = nx
            dx = L / nx

            print*, "Simulation avec nx =", nx

            allocate(U(nx+2, 3), phi_g(nx+2), P(nx+2))

            call initialize(U, dx, g, phi_g)
            call simulate(U, P, nx, g, L, tmax, phi_g)

            allocate(U_interp(nx+2, 4))
            call interpolate_solution(U_ref, U_interp, nx_ref, nx)

            allocate(Us(nx+2, 3))
            call initialize(Us, dx, g, phi_g)
            ! Construire le nom du fichier
            write(filename, '(A,I0,A)') "solution_", nx, ".txt"
            call output_results(U, Us, P, nx, L, filename)
            !call output_results(U_interp, Us, P, nx, L, filename)

            !call calculate_infinity_norm(U, U_interp, Us, nx, erreur_P, erreur_u)
            call calculate_L2_norm(U, U_interp, Us, nx, erreur_P, erreur_u)

            erreurs_P(i) = erreur_P
            erreurs_u(i) = erreur_u

            print*, "Erreur en norme infinie pour nx =", nx
            print*, "Erreur P :", erreur_P, ", Erreur u :", erreur_u

            deallocate(U_interp, Us)
            deallocate(U, phi_g, P)
        end do

        open(unit=20, file="erreurs.txt", status="unknown", action="write")
        write(20, '(A)') "nx      Erreur_P    Erreur_u"
        do i = 1, 5
            write(20, '(I4, F12.6, F12.6)') nx_values(i), erreurs_P(i), erreurs_u(i)
        end do
        close(20)

        print*, "Test de convergence terminé. Résultats enregistrés dans 'erreurs.txt'."
    end if

contains

    !################# Interpolation de la solution de référence pour le calcul des erreurs #################
    subroutine interpolate_solution(U_ref, U_interp, nx_ref, nx)
        real(PR), dimension(:,:), intent(in) :: U_ref
        real(PR), dimension(:,:), intent(out) :: U_interp
        integer, intent(in) :: nx_ref, nx
        integer :: i, j
        real(PR) :: ratio

        ratio = real(nx_ref) / real(nx)

        do i = 1, nx
            j = min(int((i-1) * ratio + int(ratio/2._PR)), nx_ref)
            !print*, j, "x = ", (i - 0.5) * 2._PR / nx, " x = ", ((j - 0.5) * 2._PR / nx_ref + (j + 0.5) * 2._PR / nx_ref)/2._PR
            U_interp(i+1, :) = (U_ref(j, :) + U_ref(j+1, :))/2._PR
        end do
    end subroutine interpolate_solution

    !################# Réalisation d'une simulation complète #################
    subroutine simulate(U, P, nx, g, L, tmax, phi_g)
        real(PR), dimension(:,:), intent(inout) :: U
        real(PR), dimension(:), intent(out) :: P
        real(PR), dimension(:), intent(in) :: phi_g
        integer, intent(in) :: nx
        real(PR), intent(in) :: g, L, tmax

        real(PR) :: dx, dt, t, lambda1, lambda2, lambda3, vitesse_u, c
        real(PR), dimension(:), allocatable :: max_vp
        real(PR), dimension(:,:), allocatable :: F, U_np1, F_arete, m_div_F, source
        integer :: i, j

        dx = L / nx
        allocate(max_vp(nx+2), F(nx+2, 3), U_np1(nx+2, 3), F_arete(nx+1, 3), m_div_F(nx, 3), source(nx, 3))

        ! Initialisation du temps
        t = 0.0
        j = 0

        ! Boucle temporelle
        do while (t < tmax)

            ! Ajout d'une perturbation sinusoïdale
            U(2, 2) = U(2, 1) * (0.1_PR * sin(8.0_PR * pi * t))

            ! Calcul de la pression
            P(2:nx+1) = (gamma - 1._PR) * (U(2:nx+1, 3) - 0.5 * U(2:nx+1, 2)**2 / U(2:nx+1, 1))

            ! Mise à jour des cellules fantômes
            call set_ghost_cells(U, P, nx)

            ! Calcul des valeurs propres
            do i = 1, nx + 2
                c = sqrt(gamma * P(i) / U(i, 1))
                vitesse_u = U(i, 2) / U(i, 1)
                lambda1 = vitesse_u - c
                lambda2 = vitesse_u
                lambda3 = vitesse_u + c
                max_vp(i) = max(abs(lambda1), abs(lambda2), abs(lambda3))
            end do

            ! Calcul du pas de temps basé sur le critère CFL
            dt = 0.9 * dx / (2 * maxval(max_vp(:)))
            if (dt <= 0.0 .or. dt > 1.0) then
                print *, "Erreur: Pas de temps non valide, dt =", dt
                stop
            end if

            ! Calcul des flux avec le schéma de Rusanov
            call calcul_flux_maille(U, P, F)
            call calcul_flux_arete(U, F, F_arete, max_vp)
            call div_flux(F_arete, m_div_F, dx)

            ! Calcul du terme source (gravité)
            call compute_source(U, source, g)

            ! Mise à jour des variables conservées
            do i = 1, nx
                if (i == 1) then
                    U_np1(i+1, :) = U(i+1, :) + dt * m_div_F(i, :) + dt * 0.5 * (source(i, :) * (phi_g(i+2) - phi_g(i+1)) / dx + source(i, :) * (phi_g(i+1) - phi_g(i)) / dx)
                else
                    U_np1(i+1, :) = U(i+1, :) + dt * m_div_F(i, :) + dt * 0.5 * (source(i, :) * (phi_g(i+2) - phi_g(i+1)) / dx + source(i-1, :) * (phi_g(i+1) - phi_g(i)) / dx)
                end if
            end do

            ! Mise à jour des variables
            U = U_np1
            t = t + dt
            j = j + 1
        end do

        ! Affichage du nombre d'itérations en temps
        print*, "Nombre d'itérations en temps : ", j

        ! Libération de la mémoire
        deallocate(max_vp, F, U_np1, F_arete, m_div_F, source)

    end subroutine simulate
  
end program euler_rusanov_gravity
  
  