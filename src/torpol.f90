!===============================================================
!  Module pour la décomposition toroïdal / poloidal
!  u = curl curl(U r) + curl(W r)
!===============================================================
module mod_torpol
  use mod_kinds
  use mod_params
  use mod_grid
  use mod_spectral_transforms
  implicit none

contains

  !-----------------------------------------------------------
  ! Calcul de u(r,θ,φ) à partir de U_hat,W_hat en spectral
  !-----------------------------------------------------------
  subroutine compute_velocity_from_UW(U_hat, W_hat, ur, uth, uph)
    complex(dp), intent(in)  :: U_hat(0:LMAX,-MMAX:MMAX,NR)
    complex(dp), intent(in)  :: W_hat(0:LMAX,-MMAX:MMAX,NR)
    real(dp),    intent(out) :: ur(0:LMAX,0:2*MMAX,NR)
    real(dp),    intent(out) :: uth(0:LMAX,0:2*MMAX,NR)
    real(dp),    intent(out) :: uph(0:LMAX,0:2*MMAX,NR)

    ! TODO:
    ! 1) Construire u_lm^r, u_lm^θ, u_lm^φ en spectral en appliquant
    !    les opérateurs vectoriels aux potentiels U_hat, W_hat
    !    (vector spherical harmonics).
    ! 2) Appeler to_physical_scalar pour chaque composante.
    !
    ! Pour l’instant, stub nul:
    ur  = 0.0_dp
    uth = 0.0_dp
    uph = 0.0_dp

  end subroutine compute_velocity_from_UW

  !-----------------------------------------------------------
  ! Calcul de la vorticité ω = curl u (en physique)
  !-----------------------------------------------------------
  subroutine compute_vorticity(ur, uth, uph, wr, wth, wph)
    real(dp), intent(in)  :: ur(0:LMAX,0:2*MMAX,NR)
    real(dp), intent(in)  :: uth(0:LMAX,0:2*MMAX,NR)
    real(dp), intent(in)  :: uph(0:LMAX,0:2*MMAX,NR)
    real(dp), intent(out) :: wr(0:LMAX,0:2*MMAX,NR)
    real(dp), intent(out) :: wth(0:LMAX,0:2*MMAX,NR)
    real(dp), intent(out) :: wph(0:LMAX,0:2*MMAX,NR)

    ! TODO: implémenter curl(u) en coordonnées sphériques
    wr  = 0.0_dp
    wth = 0.0_dp
    wph = 0.0_dp
  end subroutine compute_vorticity

  !-----------------------------------------------------------
  ! Terme non-linéaire NL = u × ω en physique
  !-----------------------------------------------------------
  subroutine compute_nonlinear_term(ur,uth,uph, wr,wth,wph, NLr,NLth,NLph)
    real(dp), intent(in)  :: ur(0:LMAX,0:2*MMAX,NR)
    real(dp), intent(in)  :: uth(0:LMAX,0:2*MMAX,NR)
    real(dp), intent(in)  :: uph(0:LMAX,0:2*MMAX,NR)
    real(dp), intent(in)  :: wr(0:LMAX,0:2*MMAX,NR)
    real(dp), intent(in)  :: wth(0:LMAX,0:2*MMAX,NR)
    real(dp), intent(in)  :: wph(0:LMAX,0:2*MMAX,NR)
    real(dp), intent(out) :: NLr(0:LMAX,0:2*MMAX,NR)
    real(dp), intent(out) :: NLth(0:LMAX,0:2*MMAX,NR)
    real(dp), intent(out) :: NLph(0:LMAX,0:2*MMAX,NR)
    integer :: it, ip, ir

    do ir = 1, NR
      do it = 0, LMAX
        do ip = 0, 2*MMAX
          NLr (it,ip,ir) = uth(it,ip,ir)*wph(it,ip,ir) - uph(it,ip,ir)*wth(it,ip,ir)
          NLth(it,ip,ir) = uph(it,ip,ir)*wr (it,ip,ir) - ur (it,ip,ir)*wph(it,ip,ir)
          NLph(it,ip,ir) = ur (it,ip,ir)*wth(it,ip,ir) - uth(it,ip,ir)*wr (it,ip,ir)
        enddo
      enddo
    enddo
  end subroutine compute_nonlinear_term

end module mod_torpol
