program precessing_sphere_dns
  use mod_kinds
  use mod_params
  use mod_grid
  use mod_field
  use mod_spectral_transforms
  use mod_torpol
  use mod_time_integration
  implicit none

  integer :: it
  logical :: first_step
  ! champs physiques pour NL
  real(dp), allocatable :: ur(:,:,:), uth(:,:,:), uph(:,:,:)
  real(dp), allocatable :: wr(:,:,:), wth(:,:,:), wph(:,:,:)
  real(dp), allocatable :: NLr(:,:,:), NLth(:,:,:), NLph(:,:,:)
  real(dp), allocatable :: U_phys(:,:,:), W_phys(:,:,:)

  call init_grids()
  call allocate_fields()

  allocate(ur(0:LMAX,0:2*MMAX,NR))
  allocate(uth(0:LMAX,0:2*MMAX,NR))
  allocate(uph(0:LMAX,0:2*MMAX,NR))
  allocate(wr(0:LMAX,0:2*MMAX,NR))
  allocate(wth(0:LMAX,0:2*MMAX,NR))
  allocate(wph(0:LMAX,0:2*MMAX,NR))
  allocate(NLr(0:LMAX,0:2*MMAX,NR))
  allocate(NLth(0:LMAX,0:2*MMAX,NR))
  allocate(NLph(0:LMAX,0:2*MMAX,NR))
  allocate(U_phys(0:LMAX,0:2*MMAX,NR))
  allocate(W_phys(0:LMAX,0:2*MMAX,NR))

  !-----------------------------------------------------------
  ! Initialisation : rotation solide approx sur W, U=0
  ! Ici en physique, puis projection en spectral
  !-----------------------------------------------------------
  U_phys = 0.0_dp
  W_phys = 0.0_dp

  ! Exemple simple : W(r=1,θ,φ) = sinθ cosφ ; on propage en r
  call init_W_physical(W_phys)

  call to_spectral_scalar(U_phys, U_hat)
  call to_spectral_scalar(W_phys, W_hat)

  LapU_hat = (0.0_dp,0.0_dp)

  RHS_W_hat     = (0.0_dp,0.0_dp)
  RHS_W_hat_old = (0.0_dp,0.0_dp)
  RHS_LapU_hat     = (0.0_dp,0.0_dp)
  RHS_LapU_hat_old = (0.0_dp,0.0_dp)

  first_step = .true.

  !-----------------------------------------------------------
  ! Boucle en temps
  !-----------------------------------------------------------
  do it = 1, nsteps

    ! U_hat,W_hat -> vitesse u en physique
    call compute_velocity_from_UW(U_hat, W_hat, ur, uth, uph)
    call compute_vorticity(ur, uth, uph, wr, wth, wph)
    call compute_nonlinear_term(ur, uth, uph, wr, wth, wph, NLr,NLth,NLph)

    ! Projeter les termes non-linéaires vers RHS_W_hat, RHS_LapU_hat
    call build_RHS_from_NL(NLr,NLth,NLph, RHS_W_hat, RHS_LapU_hat)

    call time_step(first_step)
    first_step = .false.

    if (mod(it,nout)==0) then
      write(*,*) 'it = ', it
      ! TODO: sortir un champ pour visualisation, par ex u_r ou streamfunction
    endif

  enddo

contains

  subroutine init_W_physical(W_phys)
    real(dp), intent(out) :: W_phys(0:LMAX,0:2*MMAX,NR)
    integer :: ir,it,ip
    real(dp) :: rr,th_loc,ph_loc

    W_phys = 0.0_dp
    do ir = 1, NR
      rr = r(ir)
      do it = 0, LMAX
        th_loc = theta(it)
        do ip = 0, 2*MMAX
          ph_loc = phi(ip)
          ! BC type Kida-Nakayama : ~sinθ cosφ à la paroi
          W_phys(it,ip,ir) = (rr/rmax) * sin(th_loc)*cos(ph_loc)
        enddo
      enddo
    enddo
  end subroutine init_W_physical

  subroutine build_RHS_from_NL(NLr,NLth,NLph, RHS_W_hat, RHS_LapU_hat)
    real(dp),    intent(in)  :: NLr(0:LMAX,0:2*MMAX,NR)
    real(dp),    intent(in)  :: NLth(0:LMAX,0:2*MMAX,NR)
    real(dp),    intent(in)  :: NLph(0:LMAX,0:2*MMAX,NR)
    complex(dp), intent(out) :: RHS_W_hat(0:LMAX,-MMAX:MMAX,NR)
    complex(dp), intent(out) :: RHS_LapU_hat(0:LMAX,-MMAX:MMAX,NR)

    ! TODO:
    ! Projeter les contributions non-linéaires dans l’équation de W et de LapU:
    !    RHS_W ~ opérateur[Toro] appliqué à NL
    !    RHS_LapU ~ opérateur[Pol] appliqué à NL
    !
    ! Pour l’instant, on met tout à zéro
    RHS_W_hat    = (0.0_dp,0.0_dp)
    RHS_LapU_hat = (0.0_dp,0.0_dp)
  end subroutine build_RHS_from_NL

end program precessing_sphere_dns
