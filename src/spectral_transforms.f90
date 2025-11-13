!===============================================================
!  Module pour les transformées sphériques / physiques
!  (placeholders pour SPHEREPACK/SHTns/FFTW)
!===============================================================
module mod_spectral_transforms
  use mod_kinds
  use mod_params
  use mod_grid
  implicit none

contains

  !-----------------------------------------------------------
  ! Transformée scalaires : f(θ,φ,r) → f_lm(r)
  ! Ici: placeholders; à remplacer par appels SPHEREPACK ou SHTns.
  !-----------------------------------------------------------
  subroutine to_spectral_scalar(f_phys, f_lm)
    ! f_phys:  dimension (0:LMAX, 0:2*MMAX, NR)
    ! f_lm:    dimension (0:LMAX, -MMAX:MMAX, NR)
    real(dp),    intent(in)  :: f_phys(0:LMAX,0:2*MMAX,NR)
    complex(dp), intent(out) :: f_lm(0:LMAX,-MMAX:MMAX,NR)
    integer :: ell, m, ir

    ! TODO: remplacer par vrais appels aux routines de transfo.
    f_lm = (0.0_dp, 0.0_dp)

    ! Exemple ultra brut de projection en φ par FFT naïve (à remplacer par FFTW):
    ! Ici juste un squelette vide.
    do ir = 1, NR
      do ell = 0, LMAX
        do m = -MMAX, MMAX
          f_lm(ell,m,ir) = (0.0_dp,0.0_dp)  ! TODO
        enddo
      enddo
    enddo

  end subroutine to_spectral_scalar

  !-----------------------------------------------------------
  ! Transformée inverse : f_lm(r) → f(θ,φ,r)
  !-----------------------------------------------------------
  subroutine to_physical_scalar(f_lm, f_phys)
    complex(dp), intent(in)  :: f_lm(0:LMAX,-MMAX:MMAX,NR)
    real(dp),    intent(out) :: f_phys(0:LMAX,0:2*MMAX,NR)
    integer :: ell, m, ir

    ! TODO: remplacer par vrais appels aux routines de transfo.
    f_phys = 0.0_dp

    do ir = 1, NR
      do ell = 0, LMAX
        do m = -MMAX, MMAX
          ! TODO : reconstruction par somme sur m,ℓ avec Y_ℓ^m(θ,φ)
        enddo
      enddo
    enddo

  end subroutine to_physical_scalar

  !-----------------------------------------------------------
  ! Version vectorielle : u_r,u_θ,u_φ <-> U,W
  ! (on pourrait passer par scalaires intermédiaires)
  !-----------------------------------------------------------

end module mod_spectral_transforms
