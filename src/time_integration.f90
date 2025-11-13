module time_integration
  use mod_kinds
  use mod_params
  use mod_grid
  use mod_poisson_radial
  use mod_operators
  use mod_torpol
  use mod_spectral_transforms
  implicit none

  ! History for Adams–Bashforth 2
  complex(dp), allocatable :: NL_W_old(:,:,:)
  complex(dp), allocatable :: NL_LapU_old(:,:,:)

contains

!============================================================
! Allocate history arrays
!============================================================
subroutine init_time_integrator()
  allocate(NL_W_old(0:LMAX,-MMAX:MMAX,NR))
  allocate(NL_LapU_old(0:LMAX,-MMAX:MMAX,NR))
  NL_W_old    = (0.0_dp,0.0_dp)
  NL_LapU_old = (0.0_dp,0.0_dp)
end subroutine init_time_integrator


!============================================================
! Build RHS of CN/AB2 equations
!============================================================
subroutine build_RHS(U_hat, W_hat, NL_W, NL_LapU, RHS_W, RHS_LapU)
  complex(dp), intent(in)  :: U_hat(0:LMAX,-MMAX:MMAX,NR)
  complex(dp), intent(in)  :: W_hat(0:LMAX,-MMAX:MMAX,NR)
  complex(dp), intent(in)  :: NL_W(0:LMAX,-MMAX:MMAX,NR)
  complex(dp), intent(in)  :: NL_LapU(0:LMAX,-MMAX:MMAX,NR)
  complex(dp), intent(out) :: RHS_W(0:LMAX,-MMAX:MMAX,NR)
  complex(dp), intent(out) :: RHS_LapU(0:LMAX,-MMAX:MMAX,NR)

  integer :: l,m,ir
  real(dp) :: alpha   ! ν dt / 2
  alpha = dt/(2.0_dp*Re)

  RHS_W     = (0.0_dp,0.0_dp)
  RHS_LapU  = (0.0_dp,0.0_dp)

  do ir = 1, NR
    do l = 0, LMAX
      do m = -MMAX, MMAX

        ! CN part: (I + α L) * W
        RHS_W(l,m,ir) = (1.0_dp,0.0_dp) * W_hat(l,m,ir)

        ! same for LapU
        RHS_LapU(l,m,ir) = (1.0_dp,0.0_dp) * ( Laplacian_of_U(l,m,ir) )

        ! AB2 nonlinear part
        RHS_W(l,m,ir) = RHS_W(l,m,ir) + &
           (1.5_dp)*NL_W(l,m,ir) - 0.5_dp*NL_W_old(l,m,ir)

        RHS_LapU(l,m,ir) = RHS_LapU(l,m,ir) + &
           (1.5_dp)*NL_LapU(l,m,ir) - 0.5_dp*NL_LapU_old(l,m,ir)

      enddo
    enddo
  enddo

end subroutine build_RHS


!============================================================
! Apply linear operators (Coriolis + geometry)
!============================================================
subroutine apply_linear_terms(W_hat, U_hat, RHS_W, RHS_LapU)
  complex(dp), intent(in)    :: W_hat(0:LMAX,-MMAX:MMAX,NR)
  complex(dp), intent(in)    :: U_hat(0:LMAX,-MMAX:MMAX,NR)
  complex(dp), intent(inout) :: RHS_W(0:LMAX,-MMAX:MMAX,NR)
  complex(dp), intent(inout) :: RHS_LapU(0:LMAX,-MMAX:MMAX,NR)

  call apply_linear_W(W_hat, U_hat, RHS_W)
  call apply_linear_LapU(LapU_hat, U_hat, W_hat, RHS_LapU)
end subroutine apply_linear_terms


!============================================================
! Time step update (CN for linear + AB2 for NL)
!============================================================
subroutine time_step(U_hat, W_hat)

  complex(dp), intent(inout) :: U_hat(0:LMAX,-MMAX:MMAX,NR)
  complex(dp), intent(inout) :: W_hat(0:LMAX,-MMAX:MMAX,NR)

  complex(dp), allocatable :: RHS_W(:,:,:)
  complex(dp), allocatable :: RHS_LapU(:,:,:)
  complex(dp), allocatable :: NL_W(:,:,:)
  complex(dp), allocatable :: NL_LapU(:,:,:)

  integer :: l,m,ir

  allocate(RHS_W(0:LMAX,-MMAX:MMAX,NR))
  allocate(RHS_LapU(0:LMAX,-MMAX:MMAX,NR))
  allocate(NL_W(0:LMAX,-MMAX:MMAX,NR))
  allocate(NL_LapU(0:LMAX,-MMAX:MMAX,NR))

  !-----------------------------------
  ! 1. Velocity & Vorticity
  !-----------------------------------
  call compute_velocity(U_hat, W_hat)
  call compute_vorticity()

  !-----------------------------------
  ! 2. NL term
  !-----------------------------------
  call compute_nonlinear()
  call project_nonlinear_to_spectral(NLr,NLth,NLph, NL_W, NL_LapU)

  !-----------------------------------
  ! 3. Build RHS
  !-----------------------------------
  call build_RHS(U_hat, W_hat, NL_W, NL_LapU, RHS_W, RHS_LapU)

  !-----------------------------------
  ! 4. Apply linear operators
  !-----------------------------------
  call apply_linear_terms(W_hat, U_hat, RHS_W, RHS_LapU)

  !-----------------------------------
  ! 5. Solve Poisson-type radial inversions
  !-----------------------------------
  do ir = 1, NR
    do l = 0, LMAX
      do m = -MMAX, MMAX
        W_hat(l,m,ir)    = solve_mode_radial_W( l, m, ir, RHS_W(l,m,ir) )
        LapU_hat(l,m,ir) = solve_mode_radial_LapU( l, m, ir, RHS_LapU(l,m,ir) )
      enddo
    enddo
  enddo

  !-----------------------------------
  ! 6. Update history for AB2
  !-----------------------------------
  NL_W_old    = NL_W
  NL_LapU_old = NL_LapU

end subroutine time_step

! Dummy wrappers for clarity
pure function solve_mode_radial_W(l,m,ir,rhs) result(sol)
  integer, intent(in) :: l,m,ir
  complex(dp), intent(in) :: rhs
  complex(dp) :: sol
  call solve_radial_W(l, (/rhs/), (/sol/) )
end function

pure function solve_mode_radial_LapU(l,m,ir,rhs) result(sol)
  integer, intent(in) :: l,m,ir
  complex(dp), intent(in) :: rhs
  complex(dp) :: sol
  call solve_radial_LapU(l, (/rhs/), (/sol/) )
end function

end module mod_time_integration
