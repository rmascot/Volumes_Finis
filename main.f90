program euler_rusanov_gravity

    use const
    use visualisation
    use functions

    implicit none
  
    ! Variables pour les calculs
    real(PR) :: dx, dt, t, L, g, tmax, lambda1, lambda2, lambda3, vitesse_u, t_sound, c
    real(PR), dimension(:), allocatable :: P, max_vp, phi_g
    !real(PR), dimension(3) :: Flux_gauche, Flux_droite
    real(PR), dimension(:,:), allocatable :: F, U, U_np1, m_div_F, U_analytical, Us
    real(PR), dimension(:,:), allocatable :: F_arete
    real(PR), dimension(:,:), allocatable :: source                ! Terme source (gravité)

    integer :: i, nx, j
  
    ! Initialisation des variables
    call read_parameters("parameters.txt", nx, g, L, tmax)

    ! allocate(U(nx, 3), source(nx, 3))
    ! allocate(U_analytical(nx, 3))
    ! allocate(P(nx), max_vp(nx), F_arete(nx-1, 3), F(nx, 3), U_np1(nx, 3), m_div_F(nx, 3))
    allocate(U(nx+2, 3), source(nx, 3), phi_g(nx+2))
    allocate(U_analytical(nx+2, 3))
    allocate(P(nx+2), max_vp(nx+2), F_arete(nx+1, 3), F(nx+2, 3), U_np1(nx+2, 3), m_div_F(nx, 3))

    print*, nx, g, L, tmax
    dx = L / nx
    call initialize(U, dx, g, phi_g)
  
    ! Initialisation du temps
    t = 0.0
    j = 0
    t_sound = 0._PR

    ! Boucle temporelle
    do while (t < tmax)
    !do j = 1, 2

        U(2,2) = U(2,1)*((1e-1)*sin(8*pi*t))

        ! Calcul de la pression
        P(2:nx+1) = (gamma - 1._PR) * (U(2:nx+1,3) - 0.5 * U(2:nx+1, 2)**2 / U(2:nx+1, 1))
        !P(:) = (gamma - 1._PR) * (U(:,3) - 0.5 * U(:, 2)**2 / U(:, 1))

        ! Mise à jour des ghost cells
        call set_ghost_cells(U, P, nx)

        ! Calcul des max des VP sur chaque maille
        do i = 1, nx + 2
            ! Calcul de la vitesse du son c
            c = sqrt(gamma * P(i) / U(i, 1))
            
            ! Vitesse du fluide
            vitesse_u = U(i, 2) / U(i, 1)
            
            ! Valeurs propres
            lambda1 = vitesse_u - c
            lambda2 = vitesse_u
            lambda3 = vitesse_u + c
            
            ! Maximum des valeurs propres
            max_vp(i) = max(abs(lambda1), abs(lambda2), abs(lambda3))
        end do

        ! Calcul du pas de temps basé sur CFL
        dt = 0.9 * dx / (2*maxval(max_vp(:)))

        if (dt <= 0.0 .or. dt > 1.0) then
            print *, "Erreur: Pas de temps non valide, dt =", dt
            stop
        end if
    
        ! Calcul des flux avec le schéma de Rusanov
        call calcul_flux_maille (U, P, F)

        call calcul_flux_arete(U, F, F_arete, max_vp)

        call div_flux(F_arete, m_div_F, dx)
    
        ! Calcul du terme source de gravité
        call compute_source(U, source, g)

        !print*, m_div_F(1024, 2)/m_div_F(1024, 1)
        print*, (source(1024, 2)*(phi_g(1024+2)-phi_g(1024+1))/dx + source(1024, 2)*(phi_g(1024+1)-phi_g(1024))/dx)
    
        ! Mise à jour des variables conservées
        do i = 1, nx
            !U_np1(i+1, :) = U(i+1, :) + dt * m_div_F(i, :) + dt * source(i, :)
            if (i == 1) then
                U_np1(i+1, :) = U(i+1, :) + dt * m_div_F(i, :) + dt * 0.5 * (source(i, :)*(phi_g(i+2)-phi_g(i+1))/dx + source(i, :)*(phi_g(i+1)-phi_g(i))/dx)
            else 
                U_np1(i+1, :) = U(i+1, :) + dt * m_div_F(i, :) + dt * 0.5 * (source(i, :)*(phi_g(i+2)-phi_g(i+1))/dx + source(i-1, :)*(phi_g(i+1)-phi_g(i))/dx)
            end if
        end do
    
        ! Mettre à jour les variables
        U = U_np1
  
        ! Incrémenter le temps
        t = t + dt

        j = j + 1

        !print*, "t = ", t

    end do

    print*, "Nombre d'itérations en temps : ", j

    ! Calcul de la solution analytique et de l'erreur
    allocate(Us(nx+2, 3))

    call initialize(Us, dx, g, phi_g)
    call compute_analytical(U_analytical, t, nx, dx, L)
  
    ! Sauvegarde des résultats dans un fichier
    call output_results(U, Us, P, nx, L)

    deallocate(U, Us, U_analytical, source)
    deallocate(P, max_vp, F, U_np1, m_div_F, F_arete)
  
  end program euler_rusanov_gravity
  
  