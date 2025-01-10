module visualisation

    use const

    implicit none

    contains

    !################# Écriture des résultats #################
    ! subroutine output_results(U, Us, P, nx, L)
    !     real(PR), dimension(:,:), intent(in) :: U, Us
    !     real(PR), dimension(:), intent(in) :: P
    !     real(PR), intent(in) :: L
    !     integer, intent(in) :: nx
    !     integer :: i
    !     real(PR) :: x
    !     character(len=20) :: filename
    
    !     filename = "results.txt"
    !     open(unit=10, file=filename, status="replace")
    
    !     write(10,*) "x, rho, rho*u, E, delta_P, delta_u, u_init, P_init"
    !     do i = 2, nx+1
    !         x = (i - 1 - 0.5) * L / nx
    !         write(10,*) x, U(i, 1), U(i, 2), U(i, 3), P(i)-Us(i,1)**gamma, (U(i, 2)/U(i, 1)), (Us(i, 2)/Us(i, 1)), Us(i,1)**gamma
    !         !write(10,*) x, Us(i, 1), Us(i, 2)/Us(i, 1), U(i, 3), Us(i,1)**gamma, (Us(i, 2)/Us(i, 1))
    !     end do
    
    !     close(10)
    !     print *, "Résultats enregistrés dans : ", filename
    ! end subroutine output_results

    subroutine output_results(U, Us, P, nx, L, filename)
        ! Arguments d'entrée
        real(PR), dimension(:,:), intent(in) :: U, Us
        real(PR), dimension(:), intent(in) :: P
        real(PR), intent(in) :: L
        integer, intent(in) :: nx
        character(len=*), intent(in) :: filename
    
        ! Variables locales
        integer :: i
        real(PR) :: x
    
        ! Ouverture du fichier
        open(unit=10, file=filename, status="replace")
    
        ! Écriture de l'en-tête
        !write(10,*) "x, rho, rho*u, E, delta_P, delta_u, u_init, P_init"
    
        ! Boucle sur les mailles
        do i = 2, nx+1
            x = (i - 1 - 0.5) * L / nx
            !write(10,*) x, U(i, 1), U(i, 2), U(i, 3), P(i)-Us(i,1)**gamma, (U(i, 2)/U(i, 1)), (Us(i, 2)/Us(i, 1)), Us(i,1)**gamma
            write(10,*) x, P(i)-Us(i,1)**gamma, (U(i, 2)/U(i, 1))
        end do
    
        ! Fermeture du fichier
        close(10)
    
        ! Message de confirmation
        print *, "Résultats enregistrés dans : ", filename
    end subroutine output_results

end module visualisation