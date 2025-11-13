!===============================================================
!  Module pour les transformées sphériques / physiques
!  (placeholders pour SPHEREPACK/SHTns/FFTW)
!===============================================================
module mod_spectral_transforms
  use mod_kinds
  use mod_params
  use mod_grid
  implicit none

  ! SHTns library pointer
  integer(c_int) :: shtns_handle

  ! Number of spectral coefficients
  integer :: nlat, nphi, nsp

  ! Work arrays mapped to SHTns layouts
  real(dp), allocatable :: phys_tmp(:,:)    ! (nlat, nphi)
  complex(dp), allocatable :: spec_tmp(:)   ! SHTns packed lm-coeffs

contains

!===========================================================
! Initialize SHTns plan
!===========================================================
subroutine init_sht()
  use iso_c_binding
  implicit none

  interface
    function shtns_init(nlat, nphi, lmax, mmax) bind(C,name="shtns_init")
      use iso_c_binding
      integer(c_int), value :: nlat, nphi, lmax, mmax
      integer(c_int) :: shtns_init
    end function
  end interface

  nlat = LMAX+1
  nphi = 2*MMAX+1

  ! Number of lm modes in SHTns packed layout:
  nsp = (LMAX+1) * (MMAX+1)

  ! Create work arrays
  allocate( phys_tmp(nlat,nphi) )
  allocate( spec_tmp(nsp) )

  ! Call SHTns initialization
  shtns_handle = shtns_init(nlat, nphi, LMAX, MMAX)

  if (shtns_handle == 0) then
    write(*,*) "ERROR: Failed to initialize SHTns"
    stop
  endif
end subroutine init_sht


!===========================================================
! Helper: map SHTns index (l,m) to packed index
!===========================================================
pure integer function lm_index(l,m)
  integer, intent(in) :: l,m
  lm_index = m*(LMAX+1) + l + 1    ! +1 for Fortran 1-based
end function lm_index


!===========================================================
! Convert f_phys(theta,phi,ir) → spectral (l,m)
! For each radius ir.
!===========================================================
subroutine to_spectral_scalar(f_phys, f_lm)
  use iso_c_binding
  implicit none

  real(dp),    intent(in)  :: f_phys(0:LMAX,0:2*MMAX,NR)
  complex(dp), intent(inout) :: f_lm(0:LMAX,-MMAX:MMAX,NR)

  integer :: ir,l,m,idx

  interface
    subroutine sht_forward(phys, spec) bind(C, name="sht_forward")
      use iso_c_binding
      real(c_double), intent(in)  :: phys(*)
      real(c_double), intent(out) :: spec(*)
    end subroutine
  end interface

  do ir = 1, NR

    ! Copy slice into phys_tmp
    phys_tmp(:,:) = transpose( f_phys(:,:,ir) )

    ! Call SHTns forward transform
    call sht_forward( phys_tmp, spec_tmp )

    ! Map packed spec_tmp to (l,m)
    do m = 0, MMAX
      do l = 0, LMAX
        idx = lm_index(l,m)
        f_lm(l, m, ir) = cmplx( spec_tmp(idx), 0.0_dp )
      enddo
    enddo

    ! Negative m (complex conjugate symmetry)
    do m = -MMAX, -1
      do l = 0, LMAX
        f_lm(l,m,ir) = conjg( f_lm(l,-m,ir) ) * (-1.0_dp)**m
      enddo
    enddo

  enddo

end subroutine to_spectral_scalar



!===========================================================
! Convert spectral f_lm → f_phys(theta,phi)
!===========================================================
subroutine to_physical_scalar(f_lm, f_phys)
  use iso_c_binding
  implicit none

  complex(dp), intent(in)  :: f_lm(0:LMAX,-MMAX:MMAX,NR)
  real(dp),    intent(out) :: f_phys(0:LMAX,0:2*MMAX,NR)

  integer :: ir,l,m,idx

  interface
    subroutine sht_inverse(spec, phys) bind(C, name="sht_inverse")
      use iso_c_binding
      real(c_double), intent(in)  :: spec(*)
      real(c_double), intent(out) :: phys(*)
    end subroutine
  end interface

  do ir = 1, NR

    ! Pack (l,m) → spec_tmp
    spec_tmp(:) = 0.0_dp

    do m = 0, MMAX
      do l = 0, LMAX
        idx = lm_index(l,m)
        spec_tmp(idx) = real( f_lm(l,m,ir), dp )
      enddo
    enddo

    ! SHTns is real-valued spherical transform by default
    call sht_inverse( spec_tmp, phys_tmp )

    ! Copy back
    f_phys(:,:,ir) = transpose( phys_tmp(:,:) )
  enddo

end subroutine to_physical_scalar

end module mod_spectral_transforms
