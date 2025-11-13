!===============================================================
!  Module pour opérateurs linéaires (diffusion, Coriolis, couplage)
!  dans l’espace spectral
!===============================================================
module mod_operators
  use mod_kinds
  use mod_params
  use mod_grid
  implicit none

contains

  !-----------------------------------------------------------
  ! Appliquer Laplacien spectral (angulaire + radial) à un champ f̂
  ! Ici: on ne traite que la partie angulaire -ℓ(ℓ+1)/r^2,
  ! la partie radiale est discrétisée via une matrice (Poisson radial)
  !-----------------------------------------------------------
  subroutine apply_angular_laplacian(f_hat, lap_hat)
    complex(dp), intent(in)  :: f_hat(0:LMAX,-MMAX:MMAX,NR)
    complex(dp), intent(out) :: lap_hat(0:LMAX,-MMAX:MMAX,NR)
    integer :: ell, m, ir
    real(dp) :: lambda

    do ir = 1, NR
      do ell = 0, LMAX
        lambda = - real(ell*(ell+1),dp) / (r(ir)**2 + 1.0e-12_dp)
        do m = -MMAX, MMAX
          lap_hat(ell,m,ir) = lambda * f_hat(ell,m,ir)
        enddo
      enddo
    enddo
  end subroutine apply_angular_laplacian

  !-----------------------------------------------------------
  ! Terme de Coriolis dans l’équation de Ŵ :
  !  -2 ε ∂W/∂φ → -2 ε (i m) Ŵ_{ℓ m}
  !-----------------------------------------------------------
  subroutine add_coriolis_W(W_hat, RHS_W_hat)
    complex(dp), intent(in)    :: W_hat(0:LMAX,-MMAX:MMAX,NR)
    complex(dp), intent(inout) :: RHS_W_hat(0:LMAX,-MMAX:MMAX,NR)
    integer :: ell, m, ir

    do ir = 1, NR
      do ell = 0, LMAX
        do m = -MMAX, MMAX
          RHS_W_hat(ell,m,ir) = RHS_W_hat(ell,m,ir) - (0.0_dp,2.0_dp*eps*real(m,dp)) * W_hat(ell,m,ir)
        enddo
      enddo
    enddo
  end subroutine add_coriolis_W

  !-----------------------------------------------------------
  ! Idem pour ∇²U (dans l’équation poloidale)
  !  -2 ε ∂(∇²U)/∂φ → -2 ε (i m) LapÛ_{ℓ m}
  !-----------------------------------------------------------
  subroutine add_coriolis_LapU(LapU_hat, RHS_LapU_hat)
    complex(dp), intent(in)    :: LapU_hat(0:LMAX,-MMAX:MMAX,NR)
    complex(dp), intent(inout) :: RHS_LapU_hat(0:LMAX,-MMAX:MMAX,NR)
    integer :: ell, m, ir

    do ir = 1, NR
      do ell = 0, LMAX
        do m = -MMAX, MMAX
          RHS_LapU_hat(ell,m,ir) = RHS_LapU_hat(ell,m,ir) - (0.0_dp,2.0_dp*eps*real(m,dp)) * LapU_hat(ell,m,ir)
        enddo
      enddo
    enddo
  end subroutine add_coriolis_LapU

  !-----------------------------------------------------------
  ! TODO: ajouter ici les termes de couplage géométrique entre U et W
  ! (ceux avec sinθ, cosθ, ∂θ dans les équations (9) et (11)),
  ! sous forme d’opérateurs linéaires en spectral.
  !-----------------------------------------------------------

end module mod_operators
