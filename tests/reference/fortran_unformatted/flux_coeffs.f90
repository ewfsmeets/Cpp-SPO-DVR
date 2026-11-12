! Writes flux_coeffs.bin: gfortran -fdefault-real-8 flux_coeffs.f90 && ./a.out
! A flux coefficient file as FluxAnalyzor.f90 writeFluxCoeffs writes it, for 3 blocks (nrotw*nx*ny),
! 2 planes and 4 energies: one record phi(:, iz, ien) per (iz, ien), ien outer. The value of
! (iredist, iz, ien) is (v, -1/v) with v = iredist + 10*iz + 100*ien (1-based).
program flux_coeffs
  implicit none
  integer, parameter :: funit = 70, nredist = 3, ncut = 2, nen = 4
  integer :: iredist, iz, ien
  real :: v
  complex :: phi(nredist, ncut, nen)
  do ien = 1, nen
    do iz = 1, ncut
      do iredist = 1, nredist
        v = REAL(iredist + 10*iz + 100*ien)
        phi(iredist, iz, ien) = CMPLX(v, -1.0/v)
      enddo
    enddo
  enddo
  open(unit=funit, file='flux_coeffs.bin', form='unformatted', status='replace', action='write')
  do ien = 1, nen
    do iz = 1, ncut
      write(funit) phi(:, iz, ien)
    enddo
  enddo
  close(funit)
end program flux_coeffs
