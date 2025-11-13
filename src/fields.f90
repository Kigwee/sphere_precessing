!===============================================================
!  Structures de données spectrales pour U et W
!===============================================================
module mod_fields
  use mod_kinds
  use mod_params
  implicit none

  ! Indexation: (ell, m, ir)
  ! ell = 0..LMAX
  ! m   = -ell..ell  (mais on stocke de -MMAX..MMAX en pratique)
  ! ir  = 1..NR

  complex(dp), allocatable :: U_hat(:,:,:)
  complex(dp), allocatable :: W_hat(:,:,:)
  complex(dp), allocatable :: U_hat_old(:,:,:)
  complex(dp), allocatable :: W_hat_old(:,:,:)

  ! Pour ∇² U (double Poisson)
  complex(dp), allocatable :: LapU_hat(:,:,:)
  complex(dp), allocatable :: LapU_hat_old(:,:,:)

  ! Terme RHS non-linéaire + Coriolis + couplage
  complex(dp), allocatable :: RHS_W_hat(:,:,:)
  complex(dp), allocatable :: RHS_W_hat_old(:,:,:)
  complex(dp), allocatable :: RHS_LapU_hat(:,:,:)
  complex(dp), allocatable :: RHS_LapU_hat_old(:,:,:)

contains

  subroutine allocate_fields()
    integer :: n_ell, n_m
    n_ell = LMAX+1
    n_m   = 2*MMAX+1

    allocate(U_hat  (n_ell, n_m, NR))
    allocate(W_hat  (n_ell, n_m, NR))
    allocate(U_hat_old  (n_ell, n_m, NR))
    allocate(W_hat_old  (n_ell, n_m, NR))

    allocate(LapU_hat (n_ell, n_m, NR))
    allocate(LapU_hat_old (n_ell, n_m, NR))

    allocate(RHS_W_hat (n_ell, n_m, NR))
    allocate(RHS_W_hat_old (n_ell, n_m, NR))
    allocate(RHS_LapU_hat (n_ell, n_m, NR))
    allocate(RHS_LapU_hat_old (n_ell, n_m, NR))

    U_hat         = (0.0_dp,0.0_dp)
    W_hat         = (0.0_dp,0.0_dp)
    U_hat_old     = (0.0_dp,0.0_dp)
    W_hat_old     = (0.0_dp,0.0_dp)
    LapU_hat      = (0.0_dp,0.0_dp)
    LapU_hat_old  = (0.0_dp,0.0_dp)
    RHS_W_hat     = (0.0_dp,0.0_dp)
    RHS_W_hat_old = (0.0_dp,0.0_dp)
    RHS_LapU_hat      = (0.0_dp,0.0_dp)
    RHS_LapU_hat_old  = (0.0_dp,0.0_dp)
  end subroutine allocate_fields

end module mod_fields
