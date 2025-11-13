!===============================================================
!  Module pour l’intégration en temps (CN + AB2)
!===============================================================
module mod_time_integration
  use mod_kinds
  use mod_params
  use mod_fields
  use mod_operators
  use mod_poisson_radial
  implicit none

contains

  subroutine time_step(first_step)
    logical, intent(in) :: first_step
    integer :: ell, m, ir
    complex(dp) :: lapW(0:LMAX,-MMAX:MMAX,NR)
    complex(dp) :: lapLapU(0:LMAX,-MMAX:MMAX,NR)
    complex(dp) :: RHS_loc(NR)
    complex(dp) :: sol(NR)
    real(dp) :: factor_CN

    call apply_angular_laplacian(W_hat, lapW)
    ! pour LapU, on aurait aussi lapLapU = Lap( LapU_hat ), si on traite 4e ordre complet

    ! Construire les RHS pour W_hat et LapU_hat (non-linéaire + Coriolis + couplage)
    ! Ici on suppose que RHS_W_hat et RHS_LapU_hat ont déjà reçu les contributions NL et couplage.
    call add_coriolis_W(W_hat, RHS_W_hat)
    call add_coriolis_LapU(LapU_hat, RHS_LapU_hat)

    factor_CN = dt / (2.0_dp*Re)

    ! Boucle sur (ell,m), résolution radiale
    do ell = 0, LMAX
      do m = -MMAX, MMAX

        ! W: (I - factor_CN*Lap) W^{n+1} = (I + factor_CN*Lap) W^n + dt * RHS_AB2
        ! AB2:
        if (first_step) then
          ! Euler pour initialiser
          do ir = 1, NR
            RHS_loc(ir) = W_hat(ell,m,ir) + dt * RHS_W_hat(ell,m,ir)   &
                          + factor_CN * lapW(ell,m,ir)
          enddo
        else
          do ir = 1, NR
            RHS_loc(ir) = W_hat(ell,m,ir) + dt * ( 1.5_dp*RHS_W_hat(ell,m,ir) &
                                     - 0.5_dp*RHS_W_hat_old(ell,m,ir) ) &
                          + factor_CN * lapW(ell,m,ir)
          enddo
        endif

        call solve_radial_W(ell, RHS_loc, sol)
        do ir = 1, NR
          W_hat(ell,m,ir) = sol(ir)
        enddo

        ! Idem pour LapU
        if (first_step) then
          do ir = 1, NR
            RHS_loc(ir) = LapU_hat(ell,m,ir) + dt * RHS_LapU_hat(ell,m,ir)
          enddo
        else
          do ir = 1, NR
            RHS_loc(ir) = LapU_hat(ell,m,ir) + dt * ( 1.5_dp*RHS_LapU_hat(ell,m,ir) &
                                        - 0.5_dp*RHS_LapU_hat_old(ell,m,ir) )
          enddo
        endif

        call solve_radial_LapU(ell, RHS_loc, sol)
        do ir = 1, NR
          LapU_hat(ell,m,ir) = sol(ir)
        enddo

      enddo
    enddo

    ! Double Poisson: ∇² U = LapU_hat → résoudre pour U_hat
    do ell = 0, LMAX
      do m = -MMAX, MMAX
        do ir = 1, NR
          ! TODO: remplacer par solve_radial_Poisson(ell, LapU_hat(ell,m,:), U_hat(ell,m,:))
          U_hat(ell,m,ir) = LapU_hat(ell,m,ir)   ! stub
        enddo
      enddo
    enddo

    ! Shift des RHS pour AB2
    RHS_W_hat_old    = RHS_W_hat
    RHS_LapU_hat_old = RHS_LapU_hat

  end subroutine time_step

end module mod_time_integration
