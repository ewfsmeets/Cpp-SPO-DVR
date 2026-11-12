! Flux check 2 (docs/validation_flux_replay/README.md): the reference Fortran's flux accumulation
! (storeFluxCoeffs) and probability (computeFluxProbability) on a small synthetic wave function,
! printed as hex bits for flux_f2_check.cpp (stdout); writes the flux files in its working directory.
! Built and run by run_checks.sh.
program flux_f2_check
  use Constants
  use GlobalData
  use MyFunctions
  use WaveFunction
  use FluxAnalyzor
  use CpuTimings
  implicit none
  integer, parameter :: nsteps = 5, ncut = 2, irc = 3
  integer :: istep, irot, ixy, ir, iz, izr, ien, iredist, nredist, nen
  real :: t, x
  complex :: rec(4)
  type (FluxOperator) :: flux

  call SetupConstantsModule()
  call InitCpuTimes()
  nrotw = 2
  nx = 1
  ny = 2
  nr = 8
  delr = 0.15
  delz = 0.1
  delx = 0.8
  dely = 0.9
  skewingAngle = Constants_PI / 3.0
  pointGroup = "C6v"
  rmass = (1.0078 + 1.0078) * camau
  rdmass = (1.0078 * 1.0078) / (1.0078 + 1.0078) * camau
  EMinWP = 0.2
  EMaxWP = 0.6
  z0 = 11.44
  dtAnalysis = 10.0
  edifin = 0.0123
  eviroi = 0.0098765
  newcalc = .TRUE.
  allocate( rgrid(nr), xgrid(nx), ygrid(ny) )
  rgrid = 0.0
  xgrid = 0.0
  ygrid = 0.0
  allocate( ScatWave(nrotw, nx*ny, nr*ncut) )
  write(*,'(A,4(1X,Z16.16))') 'G', transfer(rmass, 0_8), transfer(rdmass, 0_8), transfer(skewingAngle, 0_8), &
                                    transfer(dtAnalysis, 0_8)

  call setupFluxAnalyzor( flux, 0.2, 0.6, 0.15, irc, ncut, "flux_reaction.dat", "flux_density.dat", &
                          "flux_coeffs.dat", "flux_coeffs_deriv.dat" )
  do istep = 0, nsteps - 1
    t = 10.0 * REAL(istep + 1)
    do iz = 1, ncut
      do ir = 1, nr
        izr = ir + nr * (iz - 1)
        do ixy = 1, nx*ny
          do irot = 1, nrotw
            x = 0.37*REAL(irot + 2*ixy) + 1.91*REAL(iz) + 0.73*REAL(istep) + 0.29*REAL(ir)
            ScatWave(irot, ixy, izr) = CMPLX(SIN(x), COS(1.7*x + 0.4))
            write(*,'(A,5(1X,I4),2(1X,Z16.16))') 'P', istep, irot, ixy, ir, iz, &
                 transfer(REAL(ScatWave(irot, ixy, izr)), 0_8), transfer(AIMAG(ScatWave(irot, ixy, izr)), 0_8)
          enddo
        enddo
      enddo
    enddo
    write(*,'(A,1X,I4,1X,Z16.16)') 'T', istep, transfer(t, 0_8)
    call storeFluxCoeffs( flux, t )
  enddo

  call readWriteFluxCoeffs( flux, "write" )
  close(70)
  call computeFluxProbability( flux )
  close(70)

  ! Read the coefficient files back (gfortran records: one per (iz, ien)).
  nredist = nrotw * nx * ny
  nen = 3
  open(71, file="flux_coeffs.dat", form="unformatted", status="old", action="read")
  open(72, file="flux_coeffs_deriv.dat", form="unformatted", status="old", action="read")
  do ien = 1, nen
    do iz = 1, ncut
      read(71) rec(1:nredist)
      do iredist = 1, nredist
        write(*,'(A,3(1X,I4),2(1X,Z16.16))') 'F', ien, iz, iredist, transfer(REAL(rec(iredist)), 0_8), &
             transfer(AIMAG(rec(iredist)), 0_8)
      enddo
      read(72) rec(1:nredist)
      do iredist = 1, nredist
        write(*,'(A,3(1X,I4),2(1X,Z16.16))') 'D', ien, iz, iredist, transfer(REAL(rec(iredist)), 0_8), &
             transfer(AIMAG(rec(iredist)), 0_8)
      enddo
    enddo
  enddo
end program flux_f2_check
