! Flux check 1 (docs/validation_flux_replay/README.md): b(k), the D table and the flux energy count
! of the reference Fortran objects, printed as hex bit patterns for flux_f1_check.cpp (stdout).
! Fortran: MyFunctions.f90 getBkCoeffWavePacket, WaveFunction.f90 getWaveFuncDerivAtZr,
! FluxAnalyzor.f90 setupFluxAnalyzor. Built and run by run_checks.sh.
program flux_f1_check
  use Constants
  use GlobalData
  use MyFunctions
  use WaveFunction
  use FluxAnalyzor
  implicit none
  integer :: i, ir, irc, icase
  real :: emin, emax, mass, kmin, kmax, k0, w, kk
  real :: cases(3,4)
  complex :: b
  complex, pointer :: d(:)
  type (FluxOperator) :: flux

  call SetupConstantsModule()
  ! b(k) of the hardwall_small packet (EMinWP 0.2, EMaxWP 0.6 eV, H2 masses, Z0 11.44).
  mass = (1.0078 + 1.0078) * camau
  emin = 0.2 * Constants_eV2Ha
  emax = 0.6 * Constants_eV2Ha
  kmin = SQRT(2.0 * mass * emin)
  kmax = SQRT(2.0 * mass * emax)
  k0 = (kmin + kmax) / 2.0
  w = k0 - kmin
  do i = 0, 400
    kk = kmin - w + REAL(i) * (kmax - kmin + 2.0 * w) / 400.0
    b = getBkCoeffWavePacket(kk, emin, emax, mass, 11.44)
    write(*,'(A,3(1X,Z16.16))') 'B', transfer(kk, 0_8), transfer(REAL(b), 0_8), transfer(AIMAG(b), 0_8)
  enddo

  ! The D column of ir_c = 9 (1-based) on nr = 56, delr = 0.15: unit vectors through the
  ! derivative routine give ddr(ir, irc).
  nr = 56
  delr = 0.15
  nrotw = 1
  nx = 1
  ny = 1
  irc = 9
  allocate( ScatWave(1, 1, nr) )
  do ir = 1, nr
    ScatWave = (0.0, 0.0)
    ScatWave(1, 1, ir) = (1.0, 0.0)
    d => getWaveFuncDerivAtZr( irc, 1 )
    write(*,'(A,1X,I4,1X,Z16.16)') 'D', ir, transfer(REAL(d(1)), 0_8)
  enddo

  ! The energy counts (printed by setupFluxAnalyzor: "Computing flux probabilty for N energies").
  newcalc = .TRUE.
  edifin = 0.0
  eviroi = 0.0
  allocate( rgrid(nr) )
  rgrid = 0.0
  cases(:,1) = (/ 0.2, 0.6, 0.1 /)
  cases(:,2) = (/ 0.2, 0.6, 0.08 /)
  cases(:,3) = (/ 0.2, 0.6, 0.15 /)
  cases(:,4) = (/ 0.40, 0.95, 0.025 /)
  do icase = 1, 4
    call setupFluxAnalyzor( flux, cases(1,icase), cases(2,icase), cases(3,icase), irc, 1, &
                            "r", "d", "c", "cd" )
  enddo
end program flux_f1_check
