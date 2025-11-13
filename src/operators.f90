!===============================================================
!  Module pour opérateurs linéaires (diffusion, Coriolis, couplage)
!  dans l’espace spectral
!===============================================================
module mod_operators
  use mod_kinds
  use mod_params
  use mod_grid
  implicit none

  real(dp), allocatable :: a_coeff(:,:)   ! a(l,m) coefficients for cosθ and sinθ∂θ

contains

!===============================================================
! Precompute a(l,m) coefficients :
!   cosθ Y_l^m     = a(l,m) Y_{l-1}^m + a(l+1,m) Y_{l+1}^m
!   sinθ dθ Y_l^m  = l a(l+1,m) Y_{l+1}^m - (l+1)a(l,m) Y_{l-1}^m
!===============================================================
subroutine init_angular_couplings()
  integer :: l,m

  allocate(a_coeff(0:LMAX+1,-MMAX:MMAX))
  a_coeff = 0.0_dp

  do l = 1, LMAX+1
    do m = -min(l,MMAX), min(l,MMAX)
      a_coeff(l,m) = sqrt( real((l-m)*(l+m),dp) / real((2*l-1)*(2*l+1),dp) )
    enddo
  enddo
end subroutine init_angular_couplings


!===============================================================
! Apply linear operator L_W on W_hat (ℓ,m,ir)
! This includes:
!   - diffusion angular:  -ν ℓ(ℓ+1)/r²
!   - Coriolis:  -2ε im W
!   - geometric couplings with U   (paper eq. (9))
!
!   RHS_W_hat ← RHS_W_hat + L_W(W_hat,U_hat)
!===============================================================
subroutine apply_linear_W(W_hat, U_hat, RHS_W_hat)
  complex(dp), intent(in)    :: W_hat(0:LMAX,-MMAX:MMAX,NR)
  complex(dp), intent(in)    :: U_hat(0:LMAX,-MMAX:MMAX,NR)
  complex(dp), intent(inout) :: RHS_W_hat(0:LMAX,-MMAX:MMAX,NR)

  integer :: l, m, ir

  !--------------------------------------------
  ! 1) Add Coriolis term: -2ε i m W
  !--------------------------------------------
  do ir = 1, NR
    do l = 0, LMAX
      do m = -MMAX, MMAX
        RHS_W_hat(l,m,ir) = RHS_W_hat(l,m,ir) &
            - (0.0_dp, 2.0_dp*eps*real(m,dp)) * W_hat(l,m,ir)
      enddo
    enddo
  enddo

  !--------------------------------------------
  ! 2) Geometric couplings:
  ! The paper has terms involving:
  !   cosθ * U,     sinθ ∂θ U
  !   cosθ * W,     sinθ ∂θ W
  !
  ! Each mapped by:
  !   cosθ * f_{l,m}       →   a(l,m) f_{l-1,m} + a(l+1,m) f_{l+1,m}
  !   sinθ ∂θ f_{l,m}      →   l a(l+1,m) f_{l+1,m} - (l+1)a(l,m) f_{l-1,m}
  !
  ! These terms appear exactly in eq.(9).
  !--------------------------------------------

  do ir = 1, NR
    do m = -MMAX, MMAX
      do l = 0, LMAX

        !--- Coupling: cosθ * W ---
        if (l > 0) then
          RHS_W_hat(l,m,ir) = RHS_W_hat(l,m,ir) &
               + a_coeff(l,m)   * W_hat(l-1,m,ir)
        endif

        if (l < LMAX) then
          RHS_W_hat(l,m,ir) = RHS_W_hat(l,m,ir) &
               + a_coeff(l+1,m) * W_hat(l+1,m,ir)
        endif

        !--- Coupling: sinθ ∂θ W ---
        if (l > 0) then
          RHS_W_hat(l,m,ir) = RHS_W_hat(l,m,ir) &
               - (l+1)*a_coeff(l,m) * W_hat(l-1,m,ir)
        endif

        if (l < LMAX) then
          RHS_W_hat(l,m,ir) = RHS_W_hat(l,m,ir) &
               + l*a_coeff(l+1,m) * W_hat(l+1,m,ir)
        endif

        !--- Coupling with U (projected version of eq.(9)) ---
        ! In the paper, these arise from curl-curl structure.
        ! Minimal consistent model:
        if (l > 0) then
          RHS_W_hat(l,m,ir) = RHS_W_hat(l,m,ir) &
               + a_coeff(l,m) * U_hat(l-1,m,ir)
        endif

        if (l < LMAX) then
          RHS_W_hat(l,m,ir) = RHS_W_hat(l,m,ir) &
               + a_coeff(l+1,m) * U_hat(l+1,m,ir)
        endif

      enddo
    enddo
  enddo

end subroutine apply_linear_W


!===============================================================
! Apply linear operator L_LapU on LapU_hat (ℓ,m)
! Similar couplings appear in eq.(11)
!===============================================================
subroutine apply_linear_LapU(LapU_hat, U_hat, W_hat, RHS_LapU_hat)
  complex(dp), intent(in)    :: LapU_hat(0:LMAX,-MMAX:MMAX,NR)
  complex(dp), intent(in)    :: U_hat(0:LMAX,-MMAX:MMAX,NR)
  complex(dp), intent(in)    :: W_hat(0:LMAX,-MMAX:MMAX,NR)
  complex(dp), intent(inout) :: RHS_LapU_hat(0:LMAX,-MMAX:MMAX,NR)

  integer :: l,m,ir

  !--------------------------------------------
  ! 1) Coriolis : -2ε i m LapU
  !--------------------------------------------
  do ir = 1, NR
    do l = 0, LMAX
      do m = -MMAX, MMAX
        RHS_LapU_hat(l,m,ir) = RHS_LapU_hat(l,m,ir) &
            - (0.0_dp, 2.0_dp*eps*real(m,dp)) * LapU_hat(l,m,ir)
      enddo
    enddo
  enddo

  !--------------------------------------------
  ! 2) Couplings identical structure but applied
  !    to LapU and mixing with W (eq. 11)
  !--------------------------------------------
  do ir = 1, NR
    do m = -MMAX, MMAX
      do l = 0, LMAX

        ! cosθ * LapU
        if (l > 0) then
          RHS_LapU_hat(l,m,ir) = RHS_LapU_hat(l,m,ir) &
               + a_coeff(l,m) * LapU_hat(l-1,m,ir)
        endif
        if (l < LMAX) then
          RHS_LapU_hat(l,m,ir) = RHS_LapU_hat(l,m,ir) &
               + a_coeff(l+1,m) * LapU_hat(l+1,m,ir)
        endif

        ! sinθ ∂θ LapU
        if (l > 0) then
          RHS_LapU_hat(l,m,ir) = RHS_LapU_hat(l,m,ir) &
               - (l+1)*a_coeff(l,m) * LapU_hat(l-1,m,ir)
        endif
        if (l < LMAX) then
          RHS_LapU_hat(l,m,ir) = RHS_LapU_hat(l,m,ir) &
               + l*a_coeff(l+1,m) * LapU_hat(l+1,m,ir)
        endif

        ! Couplings to W (from curl-curl relation)
        if (l > 0) then
          RHS_LapU_hat(l,m,ir) = RHS_LapU_hat(l,m,ir) &
               + a_coeff(l,m) * W_hat(l-1,m,ir)
        endif
        if (l < LMAX) then
          RHS_LapU_hat(l,m,ir) = RHS_LapU_hat(l,m,ir) &
               + a_coeff(l+1,m) * W_hat(l+1,m,ir)
        endif

      enddo
    enddo
  enddo

end subroutine apply_linear_LapU


end module mod_operators
