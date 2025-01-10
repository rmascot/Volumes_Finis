module const

    implicit none
  
    ! Paramètres du problème
    integer, parameter :: PR = 8
    !integer, parameter :: nx = 400          ! Nombre de cellules
    real(PR), parameter :: gamma = 1.4      ! Rapport des chaleurs spécifiques
    !real(PR), parameter :: c = 341.0       ! Rapport des chaleurs spécifiques
    !real(PR), parameter :: g = 0.0          ! Gravité
    !real(PR), parameter :: L = 1.0           ! Longueur du domaine
    !real(PR), parameter :: tmax = 0.001
    real(PR), parameter :: pi = 4.0*atan(1.0)
    !real(PR), parameter :: dx = L / nx

end module const