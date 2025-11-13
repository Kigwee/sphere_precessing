module mod_params
  use mod_kinds
  implicit none

  !--- paramètres physiques ---
  real(dp), parameter :: Re   = 200.0_dp   ! Reynolds
  real(dp), parameter :: eps  = 0.1_dp     ! Poincare number 

  !--- géométrie de la sphère ---
  real(dp), parameter :: rmin = 0.0_dp
  real(dp), parameter :: rmax = 1.0_dp

  !--- résolution spectrale ---
  integer, parameter :: LMAX = 32       ! max degré l
  integer, parameter :: MMAX = 32       ! max m
  integer, parameter :: NR   = 32       ! points radiaux 

  !--- temps ---
  real(dp), parameter :: dt      = 1.0e-3_dp
  integer, parameter :: nsteps   = 1000
  integer, parameter :: nout     = 100

end module mod_params
