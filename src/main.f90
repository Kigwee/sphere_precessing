program dns_sphere
  use mod_kinds
  use mod_params
  use mod_grid
  use mod_spectral_transforms
  use mod_poisson_radial
  use mod_operators
  use mod_torpol
  use mod_time_integration
  implicit none

  integer :: it
  character(len=256) :: fname
  real(dp) :: t, energy

  !------------------------------------------------
  ! 0. Initialisation générale
  !------------------------------------------------
  call init_params()            ! Si tu veux mettre des valeurs ici
  call init_grid()              ! r,theta,phi
  call init_sht()               ! SHTns FFT+Legendre
  call init_radial_matrices()   ! Poisson radial CN
  call init_angular_couplings() ! a(l,m)
  call init_torpol()            ! allocate u,w arrays
  call init_time_integrator()   ! allocate NL_old

  ! Allocate main unknowns
  allocate(U_hat(0:LMAX,-MMAX:MMAX,NR))
  allocate(W_hat(0:LMAX,-MMAX:MMAX,NR))
  allocate(LapU_hat(0:LMAX,-MMAX:MMAX,NR))

  U_hat    = (0.0_dp,0.0_dp)
  W_hat    = (0.0_dp,0.0_dp)
  LapU_hat = (0.0_dp,0.0_dp)

  !------------------------------------------------
  ! 1. Conditions initiales (cas du papier)
  !    W_{1,1}(r) = sinθ cosφ = mode (l=1,m=1)
  !------------------------------------------------
  call init_conditions(U_hat, W_hat)

  t = 0.0_dp

  write(*,*) "---------------------------------------------"
  write(*,*) "   DNS tor/pol SHTns + Poisson CN/AB2 "
  write(*,*) "---------------------------------------------"
  write(*,*) "LMAX=",LMAX," MMAX=",MMAX," NR=",NR
  write(*,*) "dt=",dt," Re=",Re," eps=",eps
  write(*,*) "---------------------------------------------"

  !------------------------------------------------
  ! 2. Boucle en temps
  !------------------------------------------------
  do it = 1, nsteps
    t = t + dt

    call time_step(U_hat, W_hat)

    if (mod(it,monitor_freq)==0) then
      energy = compute_energy(U_hat,W_hat)
      write(*,'("it=",I6,"  t=",F10.4,"  E=",ES12.5)') it,t,energy
    endif

    if (mod(it,output_freq)==0) then
      write(fname,'("output/W_",I6.6,".bin")') it
      call write_field_complex(fname, W_hat)

      write(fname,'("output/U_",I6.6,".bin")') it
      call write_field_complex(fname, U_hat)
    endif
  enddo

  write(*,*) "Simulation terminée."

contains

  !----------------------------------------------------
  ! CI spécifique au papier : W = sinθ cosφ
  !----------------------------------------------------
  subroutine init_conditions(U_hat, W_hat)
    complex(dp), intent(inout) :: U_hat(0:LMAX,-MMAX:MMAX,NR)
    complex(dp), intent(inout) :: W_hat(0:LMAX,-MMAX:MMAX,NR)

    real(dp), allocatable :: W_phys(:,:,:)
    integer :: ir,i,j

    allocate(W_phys(0:LMAX,0:2*MMAX,NR))

    ! Build sinθ cosφ everywhere
    do ir = 1, NR
      do i = 0, LMAX
        do j = 0, 2*MMAX
          W_phys(i,j,ir) = sin(theta(i)) * cos(phi(j))
        enddo
      enddo
    enddo

    ! Transform into spectral
    call to_spectral_scalar(W_phys, W_hat)

    ! U initialized à 0
    U_hat = (0.0_dp,0.0_dp)
  end subroutine init_conditions


  !----------------------------------------------------
  ! Compute total kinetic energy from potentials
  ! Approximate: E = ∑ |W|² + |∇²U|²
  !----------------------------------------------------
  function compute_energy(U_hat, W_hat) result(E)
    complex(dp), intent(in) :: U_hat(0:LMAX,-MMAX:MMAX,NR)
    complex(dp), intent(in) :: W_hat(0:LMAX,-MMAX:MMAX,NR)
    real(dp) :: E
    integer :: l,m,ir

    E = 0.0_dp
    do ir=1,NR
      do l=0,LMAX
        do m=-MMAX,MMAX
          E = E + abs(W_hat(l,m,ir))**2 + abs(LapU_hat(l,m,ir))**2
        enddo
      enddo
    enddo
  end function compute_energy


  !----------------------------------------------------
  ! Writing binary field
  !----------------------------------------------------
  subroutine write_field_complex(filename, field)
    character(len=*), intent(in) :: filename
    complex(dp), intent(in) :: field(0:LMAX,-MMAX:MMAX,NR)
    integer :: u

    open(newunit=u,file=filename,form="unformatted",access="stream")
    write(u) field
    close(u)
  end subroutine write_field_complex

end program dns_sphere
