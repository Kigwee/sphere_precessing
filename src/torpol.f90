!===============================================================
!  Module pour la décomposition toroïdal / poloidal
!  u = curl curl(U r) + curl(W r)
!===============================================================
module torpol
  use mod_kinds
  use mod_params
  use mod_grid
  use mod_spectral_transforms
  implicit none

  ! Physical velocity and vorticity fields
  real(dp), allocatable :: ur(:,:,:), uth(:,:,:), uph(:,:,:)
  real(dp), allocatable :: wr(:,:,:), wth(:,:,:), wph(:,:,:)

contains
!============================================================
! Allocate physical arrays
!============================================================
subroutine init_torpol()
  allocate( ur(0:LMAX,0:2*MMAX,NR) )
  allocate( uth(0:LMAX,0:2*MMAX,NR) )
  allocate( uph(0:LMAX,0:2*MMAX,NR) )
  allocate( wr(0:LMAX,0:2*MMAX,NR) )
  allocate( wth(0:LMAX,0:2*MMAX,NR) )
  allocate( wph(0:LMAX,0:2*MMAX,NR) )
end subroutine init_torpol


!============================================================
! Compute velocity from U_hat and W_hat
!============================================================
subroutine compute_velocity(U_hat, W_hat)
  complex(dp), intent(in) :: U_hat(0:LMAX,-MMAX:MMAX,NR)
  complex(dp), intent(in) :: W_hat(0:LMAX,-MMAX:MMAX,NR)

  complex(dp), allocatable :: Ur_hat(:,:), Ut_hat(:,:), Up_hat(:,:)
  real(dp), allocatable :: u_tmp(:,:)
  integer :: ir,l,m,idx

  allocate(Ur_hat(0:LMAX,-MMAX:MMAX))
  allocate(Ut_hat(0:LMAX,-MMAX:MMAX))
  allocate(Up_hat(0:LMAX,-MMAX:MMAX))
  allocate(u_tmp(0:LMAX,0:2*MMAX))

  do ir = 1, NR

    !----- Compute radial velocity -----
    do l = 0, LMAX
      do m = -MMAX, MMAX
        Ur_hat(l,m) = real(l*(l+1),dp) / (r(ir)**2) * U_hat(l,m,ir)
      enddo
    enddo

    call to_physical_scalar(Ur_hat, u_tmp)
    ur(:,:,ir) = u_tmp(:,:)

    !----- Compute tangential velocities -----
    ! uθ and uφ use SHTns derivative transforms

    ! SHTns needs packed lm-field:
    call shtns_scalar_derivatives(U_hat(:,:,ir), W_hat(:,:,ir),  &
          Ut_hat, Up_hat, r(ir))

    call to_physical_scalar(Ut_hat, u_tmp)
    uth(:,:,ir) = u_tmp(:,:)

    call to_physical_scalar(Up_hat, u_tmp)
    uph(:,:,ir) = u_tmp(:,:)

  enddo

end subroutine compute_velocity



!============================================================
! Compute vorticity = curl(u)
!============================================================
subroutine compute_vorticity()
  integer :: ir

  do ir = 1, NR
    call curl_vector(ur(:,:,ir), uth(:,:,ir), uph(:,:,ir),  &
                     wr(:,:,ir), wth(:,:,ir), wph(:,:,ir), r(ir))
  enddo

end subroutine compute_vorticity



!============================================================
! Non-linear term: NL = u × ω
!============================================================
subroutine compute_nonlinear()
  integer :: ir,i,j

  do ir = 1, NR
    do i = 0, LMAX
      do j = 0, 2*MMAX

        NLr(i,j,ir)  = uth(i,j,ir)*wph(i,j,ir) - uph(i,j,ir)*wth(i,j,ir)
        NLth(i,j,ir) = uph(i,j,ir)*wr(i,j,ir)  - ur(i,j,ir)*wph(i,j,ir)
        NLph(i,j,ir) = ur(i,j,ir)*wth(i,j,ir)  - uth(i,j,ir)*wr(i,j,ir)

      enddo
    enddo
  enddo

end subroutine compute_nonlinear



!============================================================
! Project NL -> RHS_W and RHS_LapU (spectral)
!============================================================
subroutine project_nonlinear_to_spectral(NLr,NLth,NLph, RHS_W, RHS_LapU)
  real(dp), intent(in) :: NLr(0:LMAX,0:2*MMAX,NR)
  real(dp), intent(in) :: NLth(0:LMAX,0:2*MMAX,NR)
  real(dp), intent(in) :: NLph(0:LMAX,0:2*MMAX,NR)
  complex(dp), intent(out) :: RHS_W(0:LMAX,-MMAX:MMAX,NR)
  complex(dp), intent(out) :: RHS_LapU(0:LMAX,-MMAX:MMAX,NR)

  integer :: ir

  do ir = 1, NR

    ! Compute NL scalar = divergence-free projection of curl NL
    call nonlinear_projection_sht(NLr(:,:,ir), NLth(:,:,ir), NLph(:,:,ir), &
                                  RHS_W(:,:,ir), RHS_LapU(:,:,ir), r(ir))

  enddo

end subroutine project_nonlinear_to_spectral


subroutine shtns_scalar_derivatives(Ulm, Wlm, Utlm, Uplm, r)
  ! Compute tangential velocities (spectral) using SHTns operators
  !
  ! u_theta,lm = (1/r) * d/dr (r U_lm) * dθY_lm  + (im / (r sinθ)) W_lm Y_lm
  ! u_phi,lm   = (im/(r sinθ)) d/dr(r U_lm) Y_lm - (1/r) W_lm dθY_lm
  !
  use kinds_mod
  use params_mod
  use spectral_transforms_mod
  implicit none

  complex(dp), intent(in)  :: Ulm(0:LMAX,-MMAX:MMAX)
  complex(dp), intent(in)  :: Wlm(0:LMAX,-MMAX:MMAX)
  complex(dp), intent(out) :: Utlm(0:LMAX,-MMAX:MMAX)
  complex(dp), intent(out) :: Uplm(0:LMAX,-MMAX:MMAX)
  real(dp), intent(in)     :: r

  interface
    subroutine sht_dtheta(spec, dspec) bind(C,name="sht_dtheta")
      use iso_c_binding
      real(c_double), intent(in) :: spec(*)
      real(c_double), intent(out):: dspec(*)
    end subroutine

    subroutine sht_dphi(spec, dspec) bind(C,name="sht_dphi")
      use iso_c_binding
      real(c_double), intent(in) :: spec(*)
      real(c_double), intent(out):: dspec(*)
    end subroutine
  end interface

  real(dp), allocatable :: dth(:), dph(:)
  real(dp), allocatable :: U_r(:)
  integer :: l,m,idx

  allocate( dth(nsp), dph(nsp), U_r(nsp) )
  allocate( Utlm, Uplm )

  ! Pack Ulm into spec_tmp (already available)
  do m = 0, MMAX
    do l = 0, LMAX
      idx = lm_index(l,m)
      spec_tmp(idx) = real( Ulm(l,m) ,dp)
    enddo
  enddo

  ! dθ of U
  call sht_dtheta(spec_tmp, dth)

  ! dφ of U  (gives i m U)
  call sht_dphi(spec_tmp, dph)

  ! Compute radial derivative term: d/dr (r U)
  do m = -MMAX,MMAX
    do l = 0, LMAX
      idx = lm_index(l,abs(m))
      U_r(idx) = (r * real(Ulm(l,m),dp))   ! real only (energy in U)
    enddo
  enddo

  ! Fill Utlm and Uplm
  do m = -MMAX, MMAX
    do l = 0, LMAX

      idx = lm_index(l,abs(m))

      ! u_theta(l,m)
      Utlm(l,m) = cmplx( ( U_r(idx) / r ) * dth(idx)  &
                     + ( real(m,dp)/(r) ) * real(Wlm(l,m),dp), 0.0_dp )

      ! u_phi(l,m)
      Uplm(l,m) = cmplx( ( real(m,dp)/r ) * U_r(idx)    &
                     - (1.0_dp/r) * dth(idx), 0.0_dp )
    enddo
  enddo

  deallocate(dth,dph,U_r)

end subroutine shtns_scalar_derivatives

subroutine curl_vector(ur,uth,uph, wr,wth,wph, r)
  use kinds_mod
  use params_mod
  implicit none

  real(dp), intent(in)  :: ur(0:LMAX,0:2*MMAX)
  real(dp), intent(in)  :: uth(0:LMAX,0:2*MMAX)
  real(dp), intent(in)  :: uph(0:LMAX,0:2*MMAX)
  real(dp), intent(out) :: wr(0:LMAX,0:2*MMAX)
  real(dp), intent(out) :: wth(0:LMAX,0:2*MMAX)
  real(dp), intent(out) :: wph(0:LMAX,0:2*MMAX)
  real(dp), intent(in)  :: r

  integer :: i,j
  real(dp) :: dth_u_phi, dph_u_th, dph_u_r
  real(dp) :: sinth, costh
  real(dp) :: dth_u_r, dth_u_th, dth_u_ph

  do i = 0, LMAX
    sinth = sin(theta(i))
    do j = 0, 2*MMAX

      ! Very simple FD on θ,φ — acceptable because resolution is high

      if (i==0) then
        dth_u_phi = (uph(i+1,j)-uph(i,j))/(theta(1)-theta(0))
        dth_u_r   = (ur(i+1,j) -ur(i,j))/(theta(1)-theta(0))
      else if (i==LMAX) then
        dth_u_phi = (uph(i,j)-uph(i-1,j))/(theta(i)-theta(i-1))
        dth_u_r   = (ur(i,j) -ur(i-1,j))/(theta(i)-theta(i-1))
      else
        dth_u_phi = (uph(i+1,j)-uph(i-1,j))/(theta(i+1)-theta(i-1))
        dth_u_r   = (ur(i+1,j) -ur(i-1,j))/(theta(i+1)-theta(i-1))
      endif

      ! dφ derivative (simple centered periodic)
      if (j==0) then
        dph_u_th = (uth(i,1)-uth(i,2*MMAX))/(phi(1)-phi(0))
        dph_u_r  = (ur(i,1) -ur(i,2*MMAX))/(phi(1)-phi(0))
      else if (j==2*MMAX) then
        dph_u_th = (uth(i,0)-uth(i,j-1))/(phi(1)-phi(0))
        dph_u_r  = (ur(i,0) -ur(i,j-1))/(phi(1)-phi(0))
      else
        dph_u_th = (uth(i,j+1)-uth(i,j-1))/(phi(j+1)-phi(j-1))
        dph_u_r  = (ur(i,j+1) -ur(i,j-1))/(phi(j+1)-phi(j-1))
      endif

      ! Curl components
      wr(i,j)  = (1.0_dp/(r*sinth)) * ( dth_u_phi - dph_u_th )
      wth(i,j) = (1.0_dp/r) * ( (1.0_dp/sinth)*dph_u_r - (uph(i,j) + r*dth_u_phi)/r )
      wph(i,j) = (1.0_dp/r) * ( (uth(i,j)+r*dth_u_r)/r - dth_u_r )

    enddo
  enddo

end subroutine curl_vector

subroutine nonlinear_projection_sht(NLr,NLth,NLph, Wlm, LapUlm, r)
  use kinds_mod
  use params_mod
  use spectral_transforms_mod
  implicit none

  real(dp), intent(in)  :: NLr(0:LMAX,0:2*MMAX)
  real(dp), intent(in)  :: NLth(0:LMAX,0:2*MMAX)
  real(dp), intent(in)  :: NLph(0:LMAX,0:2*MMAX)
  complex(dp), intent(out) :: Wlm(0:LMAX,-MMAX:MMAX)
  complex(dp), intent(out) :: LapUlm(0:LMAX,-MMAX:MMAX)
  real(dp), intent(in)     :: r

  complex(dp) :: Rlm(0:LMAX,-MMAX:MMAX)
  complex(dp) :: Tlm(0:LMAX,-MMAX:MMAX)
  integer :: l,m

  ! Step 1: Transform NL_r and NL_phi to spherical harmonics
  call to_spectral_scalar(NLr,  Rlm)
  call to_spectral_scalar(NLph, Tlm)

  ! Step 2: Projection formulas from MHD tor/pol
  do m = -MMAX, MMAX
    do l = 1, LMAX
      LapUlm(l,m) = Rlm(l,m)
      Wlm(l,m)    = Tlm(l,m) / real(l*(l+1),dp)
    enddo
  enddo

end subroutine nonlinear_projection_sht

end module mod_torpol
