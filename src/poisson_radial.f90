!===============================================================
!  Poisson radial (1D en r) pour chaque (ℓ,m) via LAPACK
!===============================================================
module mod_poisson_radial
  use mod_kinds
  use mod_params
  use mod_grid
  implicit none

  ! Matrices stockées pleine (NR×NR) pour chaque ℓ
  real(dp), allocatable :: A_W(:,:,:)
  real(dp), allocatable :: A_LapU(:,:,:)

  ! Pivots LAPACK
  integer,  allocatable :: ipiv_W(:,:)
  integer,  allocatable :: ipiv_LapU(:,:)

contains

!===============================================================
! Build radial matrices for W and LapU:
!   A_W(ℓ)    = I - α L_ℓ
!   A_LapU(ℓ) = I - α L_ℓ   with modified BC at r=1
!
! L_ℓ f = f'' + (2/r) f' - ℓ(ℓ+1)/r² f
! α = dt/(2Re)
!===============================================================
subroutine init_radial_matrices()
  integer :: ell, i
  real(dp) :: dr, alpha
  real(dp) :: ri
  integer :: n

  n = NR
  dr = (rmax - rmin) / real(NR-1,dp)
  alpha = dt / (2.0_dp * Re)

  allocate(A_W(0:LMAX, NR, NR))
  allocate(A_LapU(0:LMAX, NR, NR))
  allocate(ipiv_W(0:LMAX, NR))
  allocate(ipiv_LapU(0:LMAX, NR))

  do ell = 0, LMAX

    A_W(ell,:,:) = 0.0_dp
    A_LapU(ell,:,:) = 0.0_dp

    !----------------------------------------------------------
    ! Interior FD stencil: i = 2..NR-1
    !----------------------------------------------------------
    do i = 2, NR-1
      ri = r(i)

      ! f''
      A_W(ell,i,i-1) = A_W(ell,i,i-1) + alpha * ( 1.0_dp / dr**2 - 1.0_dp/(ri*2.0_dp*dr) )
      A_W(ell,i,i)   = A_W(ell,i,i)   + alpha * ( -2.0_dp / dr**2 - ell*(ell+1)/ri**2 )
      A_W(ell,i,i+1) = A_W(ell,i,i+1) + alpha * ( 1.0_dp / dr**2 + 1.0_dp/(ri*2.0_dp*dr) )

      ! Identity part (Crank–Nicolson)
      A_W(ell,i,i) = A_W(ell,i,i) + 1.0_dp
    enddo

    !----------------------------------------------------------
    ! Boundary r=0 (regularity): W(0)=0
    !----------------------------------------------------------
    A_W(ell,1,:) = 0.0_dp
    A_W(ell,1,1) = 1.0_dp

    !----------------------------------------------------------
    ! Boundary r=1 (Dirichlet) for W: W(1)=given
    !----------------------------------------------------------
    A_W(ell,NR,:) = 0.0_dp
    A_W(ell,NR,NR) = 1.0_dp

    !----------------------------------------------------------
    ! Copy into A_LapU then change BC
    !----------------------------------------------------------
    A_LapU(ell,:,:) = A_W(ell,:,:)

    ! BC for U: U(1)=0
    A_LapU(ell,NR,:) = 0.0_dp
    A_LapU(ell,NR,NR) = 1.0_dp

    ! derivative BC: U'(1)=0 → (U_NR - U_{NR-2})/(2dr)=0
    A_LapU(ell,NR-1,:) = 0.0_dp
    A_LapU(ell,NR-1,NR)   =  1.0_dp/(2.0_dp*dr)
    A_LapU(ell,NR-1,NR-2) = -1.0_dp/(2.0_dp*dr)

  enddo

end subroutine init_radial_matrices

!===============================================================
! Solve A_W(ell)*x = rhs for complex rhs
!===============================================================
subroutine solve_radial_W(ell, rhs, sol)
  integer, intent(in) :: ell
  complex(dp), intent(in)  :: rhs(NR)
  complex(dp), intent(out) :: sol(NR)
  integer :: info, i
  real(dp) :: b_re(NR), b_im(NR)
  real(dp) :: A(NR,NR)

  ! Copy matrix
  A = A_W(ell,:,:)

  ! Real system for Re(rhs)
  do i = 1, NR
    b_re(i) = real(rhs(i),dp)
  enddo
  call dgesv(NR, 1, A, NR, ipiv_W(ell,:), b_re, NR, info)

  ! Real system for Im(rhs)
  A = A_W(ell,:,:)   ! need to reload matrix (destroyed by dgesv)
  do i = 1, NR
    b_im(i) = aimag(rhs(i))
  enddo
  call dgesv(NR, 1, A, NR, ipiv_W(ell,:), b_im, NR, info)

  do i = 1, NR
    sol(i) = cmplx(b_re(i), b_im(i), kind=dp)
  enddo
end subroutine solve_radial_W


!===============================================================
! Solve A_LapU(ell)*x = rhs
!===============================================================
subroutine solve_radial_LapU(ell, rhs, sol)
  integer, intent(in) :: ell
  complex(dp), intent(in)  :: rhs(NR)
  complex(dp), intent(out) :: sol(NR)
  integer :: info, i
  real(dp) :: b_re(NR), b_im(NR)
  real(dp) :: A(NR,NR)

  A = A_LapU(ell,:,:)

  do i = 1, NR
    b_re(i) = real(rhs(i),dp)
  enddo
  call dgesv(NR, 1, A, NR, ipiv_LapU(ell,:), b_re, NR, info)

  A = A_LapU(ell,:,:)
  do i = 1, NR
    b_im(i) = aimag(rhs(i))
  enddo
  call dgesv(NR, 1, A, NR, ipiv_LapU(ell,:), b_im, NR, info)

  do i = 1, NR
    sol(i) = cmplx(b_re(i), b_im(i), kind=dp)
  enddo
end subroutine solve_radial_LapU

end module mod_poisson_radial

