module mod_grid
  use mod_kinds
  use mod_params
  implicit none

  real(dp), allocatable :: r(:)      ! radial grid (collocation)
  real(dp), allocatable :: w_r(:)    ! radial weights (si Gauss)
  real(dp), allocatable :: theta(:)  ! grid physique 
  real(dp), allocatable :: phi(:)

contains

  subroutine init_grids()
    integer :: i
    real(dp) :: dr

    !--- radial grid simple en collocation uniforme (remplaçable par Jacobi) ---
    allocate(r(NR), w_r(NR))
    dr = (rmax - rmin) / real(NR-1,dp)
    do i = 1, NR
      r(i)   = rmin + dr * real(i-1,dp)
      w_r(i) = dr
    enddo

    !--- grilles physiques en θ,φ (optionnelles, utiles pour sortie et NL) ---
    allocate(theta(0:LMAX))   ! par exemple, LMAX+1 points en θ
    allocate(phi(0:2*MMAX))   ! 2MMAX+1 points en φ

    do i = 0, LMAX
      theta(i) = real(i,dp) * 3.14159265358979323846_dp / real(LMAX,dp)
    enddo

    do i = 0, 2*MMAX
      phi(i) = real(i,dp) * 2.0_dp*3.14159265358979323846_dp / real(2*MMAX+1,dp)
    enddo

  end subroutine init_grids

end module mod_grid
