!===============================================================
!  Poisson radial (1D en r) pour chaque (ℓ,m) via LAPACK
!===============================================================
module mod_poisson_radial
  use mod_kinds
  use mod_params
  use mod_grid
  implicit none

  ! Matrices factoriées pour chaque ℓ (et éventuellement m)
  ! Ici on simplifie: même matrice pour tous m à ℓ donné.
  real(dp), allocatable :: A_W(:,:,:)
  real(dp), allocatable :: A_LapU(:,:,:)
  integer,  allocatable :: ipiv_W(:,:), ipiv_LapU(:,:)

contains

  subroutine init_radial_matrices()
    integer :: ell, n_ell

    n_ell = LMAX+1
    allocate(A_W(n_ell, NR, NR))
    allocate(A_LapU(n_ell, NR, NR))
    allocate(ipiv_W(n_ell, NR))
    allocate(ipiv_LapU(n_ell, NR))

    ! TODO:
    ! Pour chaque ell, construire la matrice radiale discrétisant:
    !  (I - dt/(2Re) * (d²/dr² + 2/r d/dr - ell(ell+1)/r²)) pour W
    !  (I - dt/(2Re) * Laplacien 4e ordre) pour LapU
    !
    ! Puis factoriser avec dgesv/dgbsv (&co) de LAPACK.
    !
    ! Pour l’instant, on met des identités.
    A_W = 0.0_dp
    A_LapU = 0.0_dp
    do ell = 0, LMAX
      A_W(ell, :,:)    = 0.0_dp
      A_LapU(ell,:,:)  = 0.0_dp
      A_W(ell,  :, :)  = 0.0_dp
      A_LapU(ell,:, :) = 0.0_dp
      ! identité:
      A_W(ell, 1:NR,1:NR)   = 0.0_dp
      A_LapU(ell,1:NR,1:NR) = 0.0_dp
      A_W(ell, 1:NR,1:NR)   = A_W(ell,1:NR,1:NR) + diag_identity(NR)
      A_LapU(ell,1:NR,1:NR) = A_LapU(ell,1:NR,1:NR) + diag_identity(NR)
    enddo

  end subroutine init_radial_matrices

  function diag_identity(n) result(mat)
    integer, intent(in) :: n
    real(dp) :: mat(n,n)
    integer :: i
    mat = 0.0_dp
    do i = 1, n
      mat(i,i) = 1.0_dp
    enddo
  end function diag_identity

  ! Solve A_W(ell) * x = b for each (ell,m)
  subroutine solve_radial_W(ell, rhs, sol)
    integer, intent(in) :: ell
    complex(dp), intent(in)  :: rhs(NR)
    complex(dp), intent(out) :: sol(NR)

    ! TODO: utiliser dgesv sur la partie réelle/imag séparément
    sol = rhs  ! stub
  end subroutine solve_radial_W

  subroutine solve_radial_LapU(ell, rhs, sol)
    integer, intent(in) :: ell
    complex(dp), intent(in)  :: rhs(NR)
    complex(dp), intent(out) :: sol(NR)
    sol = rhs  ! stub
  end subroutine solve_radial_LapU

end module mod_poisson_radial
