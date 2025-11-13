module mod_params
  use mod_kinds
  implicit none

  !===========================================================
  ! Paramètres globaux du DNS
  !===========================================================

  ! Résolution angulaire (Harmoniques sphériques)
  integer :: LMAX      = 16
  integer :: MMAX      = 16

  ! Nombre de points radiaux
  integer :: NR        = 32

  ! Domain
  real(dp) :: rmin     = 0.0_dp
  real(dp) :: rmax     = 1.0_dp

  ! Temps
  real(dp) :: dt       = 1e-3_dp
  integer  :: nsteps   = 10000

  ! Physique
  real(dp) :: Re       = 100.0_dp     ! Reynolds
  real(dp) :: eps      = 1.0_dp       ! Coefficient de Coriolis

  ! Monitoring & output
  integer :: monitor_freq = 100
  integer :: output_freq  = 500

  !===========================================================
  ! Conteneurs principaux (ALLOC dans main)
  !===========================================================
  ! Ces trois variables sont rendues globales pour simplifier.
  ! Elles seront ALLOCATE dans main.
  complex(dp), allocatable :: U_hat(:,:,:)
  complex(dp), allocatable :: W_hat(:,:,:)
  complex(dp), allocatable :: LapU_hat(:,:,:)   ! Laplacien(U)

contains

  !===========================================================
  ! Initialisation unique
  !===========================================================
  subroutine init_params()
    implicit none

    write(*,*) "-----------------------------------------"
    write(*,*) " Initialisation des paramètres"
    write(*,*) "-----------------------------------------"
    write(*,'("  LMAX=",I4,"  MMAX=",I4,"  NR=",I4)') LMAX,MMAX,NR
    write(*,'("  dt=",ES10.3,"  nsteps=",I6)') dt,nsteps
    write(*,'("  Re=",F8.2,"  eps=",F8.2)') Re,eps
    write(*,'("  output_freq=",I6,"  monitor_freq=",I6)') output_freq, monitor_freq
    write(*,*) "-----------------------------------------"
  end subroutine init_params

end module mod_params
