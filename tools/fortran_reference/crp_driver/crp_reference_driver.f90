! Reference driver for the C++ port of the Fortran CRP library (src/crplib/).
!
! Calls the public routines of the Fortran library (objects of the reference build, see build.sh)
! on inputs given as hexadecimal IEEE bit patterns and prints the results bit for bit (hex) and
! readable (ES). The driver itself does no floating-point arithmetic and no trigonometry.
! Every run ends with the line CRP_DRIVER_DONE (a Fortran STOP inside the library exits with
! status 0, so the caller must check for it). One data set per process (the library keeps SAVEd
! module state).
!
! Usage: crp_reference_driver <mode> [arguments] < input > output
!   tables                              Fourier coefficient tables
!   wrapper   <dataset>                 points (u v Z r theta phi)        -> V (eV), parameters
!   molecule  <dataset> <length> <energy>  points (Za Xa Ya Zb Xb Yb)   -> V, F(1:6) in those units
!   atom      <dataset> <length> <energy>  points (x y z)               -> V, F(1:3)
!   gasphase  <dataset>                 points (Zcc r)                    -> Fso dFso V2D dZ dR
!   cambcord  <dataset>                 points (xa ya za xb yb zb)        -> x y z r theta phi
!   spline1d                            n ic1 icn m / cm1 cmn / x(n) / y(n) / t(m)
!   spline2d                            nx ny m / x(nx) / y(ny) / f(nx*ny) / (x y)(m)
!   spline3d                            nx ny nz m / x / y / z / f / der1(nx*ny) / der2 / (x y z)(m)
!   trispl                              n / x(n) / y(n)
!   ubicar                              n m / x(n) / t(m)
! Point input: first line the number of points, then one point per line (hex values).
! Unit codes as in UnitConversion.f90 (length 0 bohr, 1 Angstrom; energy 2 Hartree, 3 eV, ...).
PROGRAM crp_reference_driver
  USE CRP_Wrapper
  USE CRP_QD_Wrapper
  USE CRP_GasPhase
  USE UnitConversion
  USE Spline
  IMPLICIT NONE

  INTERFACE
    SUBROUTINE dump_crp_tables() BIND(C, name='dump_crp_tables')
    END SUBROUTINE
  END INTERFACE

  CHARACTER(LEN=32)  :: mode
  CHARACTER(LEN=512) :: dataset

  CALL GET_COMMAND_ARGUMENT(1, mode)
  CALL GET_COMMAND_ARGUMENT(2, dataset)

  SELECT CASE (TRIM(mode))
  CASE ('tables')
    CALL Initialize_CRP_Wrapper_Module()
    CALL dump_crp_tables()
  CASE ('wrapper')
    CALL run_wrapper()
  CASE ('molecule')
    CALL run_molecule()
  CASE ('atom')
    CALL run_atom()
  CASE ('gasphase')
    CALL run_gasphase()
  CASE ('cambcord')
    CALL run_cambcord()
  CASE ('spline1d')
    CALL run_spline1d()
  CASE ('spline2d')
    CALL run_spline2d()
  CASE ('spline3d')
    CALL run_spline3d()
  CASE ('trispl')
    CALL run_trispl()
  CASE ('ubicar')
    CALL run_ubicar()
  CASE DEFAULT
    PRINT '(A)', 'unknown mode: ' // TRIM(mode)
    STOP 1
  END SELECT
  PRINT '(A)', 'CRP_DRIVER_DONE'

CONTAINS

  ! A double from 16 hex digits (read as two 32-bit halves: no overflow of INTEGER(8)).
  REAL(8) FUNCTION from_hex(token)
    CHARACTER(LEN=*), INTENT(IN) :: token
    INTEGER(8) :: high, low
    READ(token(1:8), '(Z8)') high
    READ(token(9:16), '(Z8)') low
    from_hex = TRANSFER(IOR(ISHFT(high, 32), low), 1.0_8)
  END FUNCTION

  SUBROUTINE read_hex(values)
    REAL(8), INTENT(OUT) :: values(:)
    CHARACTER(LEN=16) :: tokens(SIZE(values))
    INTEGER :: i
    IF (SIZE(values) == 0) RETURN
    READ(*, *) tokens
    DO i = 1, SIZE(values)
      values(i) = from_hex(tokens(i))
    END DO
  END SUBROUTINE

  SUBROUTINE read_points(width, points)
    INTEGER, INTENT(IN) :: width
    REAL(8), ALLOCATABLE, INTENT(OUT) :: points(:,:)
    INTEGER :: n, i
    READ(*, *) n
    ALLOCATE(points(width, n))
    DO i = 1, n
      CALL read_hex(points(:, i))
    END DO
  END SUBROUTINE

  SUBROUTINE put(label, values)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(8), INTENT(IN) :: values(:)
    INTEGER :: i
    DO i = 1, SIZE(values)
      WRITE(*, '(A,1X,I0,1X,Z16.16,1X,ES26.17E3)') label, i, TRANSFER(values(i), 0_8), values(i)
    END DO
  END SUBROUTINE

  SUBROUTINE read_units(units_out)
    TYPE(Units), INTENT(OUT) :: units_out
    CHARACTER(LEN=16) :: text
    CALL GET_COMMAND_ARGUMENT(3, text)
    READ(text, *) units_out%LengthUnit
    CALL GET_COMMAND_ARGUMENT(4, text)
    READ(text, *) units_out%EnergyUnit
    units_out%MassUnit = UNITS_MELEC
    units_out%TimeUnit = UNITS_AUT
  END SUBROUTINE

  SUBROUTINE setup_dataset()
    CALL Set_PES_DataSet(TRIM(dataset))
    CALL Initialize_CRP_QD_Wrapper_Module()
  END SUBROUTINE

  SUBROUTINE run_wrapper()
    REAL(8) :: xlat, ylat, zasy, skew
    CHARACTER(LEN=10) :: group
    REAL(8), ALLOCATABLE :: points(:,:), v(:)
    INTEGER :: i
    CALL Set_PES_DataSet(TRIM(dataset))
    CALL Get_CRP_QD_Wrapper_Params(xlat, ylat, zasy, skew, group)
    PRINT '(A,1X,A)', 'point_group', TRIM(group)
    CALL put('params', (/ xlat, ylat, zasy, skew /))
    CALL read_points(6, points)
    ALLOCATE(v(SIZE(points, 2)))
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC)
    DO i = 1, SIZE(points, 2)
      v(i) = Get_CRP_QD_Wrapper_Potential_6D(points(1,i), points(2,i), points(3,i), points(4,i), &
                                             points(5,i), points(6,i), 0, 0)
    END DO
    !$OMP END PARALLEL DO
    DO i = 1, SIZE(points, 2)
      CALL put('V', (/ v(i) /))
    END DO
  END SUBROUTINE

  SUBROUTINE run_molecule()
    TYPE(Units) :: out_units
    REAL(8), ALLOCATABLE :: points(:,:), v(:), f(:,:)
    REAL(8) :: r(6), force(6), value
    INTEGER :: i
    CALL setup_dataset()
    CALL read_units(out_units)
    CALL read_points(6, points)
    ALLOCATE(v(SIZE(points, 2)), f(6, SIZE(points, 2)))
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(r, force, value)
    DO i = 1, SIZE(points, 2)
      r = points(:, i)
      CALL Calculate_CRP_MoleculeSurface_PES(PESList(1), r, value, force, out_units)
      v(i) = value
      f(:, i) = force
    END DO
    !$OMP END PARALLEL DO
    DO i = 1, SIZE(points, 2)
      CALL put('V', (/ v(i) /))
      CALL put('F', f(:, i))
    END DO
  END SUBROUTINE

  SUBROUTINE run_atom()
    TYPE(Units) :: out_units
    REAL(8), ALLOCATABLE :: points(:,:), v(:), f(:,:)
    REAL(8) :: r(3), force(3), value
    INTEGER :: i
    CALL setup_dataset()
    CALL read_units(out_units)
    CALL read_points(3, points)
    ALLOCATE(v(SIZE(points, 2)), f(3, SIZE(points, 2)))
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(r, force, value)
    DO i = 1, SIZE(points, 2)
      r = points(:, i)
      CALL Calculate_CRP_AtomSurface_PES(PESList(1)%AtomSurface, r, value, force, out_units)
      v(i) = value
      f(:, i) = force
    END DO
    !$OMP END PARALLEL DO
    DO i = 1, SIZE(points, 2)
      CALL put('V', (/ v(i) /))
      CALL put('F', f(:, i))
    END DO
  END SUBROUTINE

  SUBROUTINE run_gasphase()
    REAL(8), ALLOCATABLE :: points(:,:)
    REAL(8) :: fso, dfso, v2d, dz, dr
    INTEGER :: i
    CALL setup_dataset()
    CALL put('zso', (/ PESList(1)%GasPhase%ZsoMin, PESList(1)%GasPhase%ZsoMax /))
    CALL read_points(2, points)
    DO i = 1, SIZE(points, 2)
      CALL Switch_Function(points(1,i), fso, dfso, PESList(1)%GasPhase)
      CALL Calculate_CRP_GasPhase_Potential(points(1,i), points(2,i), v2d, dz, dr, PESList(1)%GasPhase)
      CALL put('G', (/ fso, dfso, v2d, dz, dr /))
    END DO
  END SUBROUTINE

  SUBROUTINE run_cambcord()
    REAL(8), ALLOCATABLE :: points(:,:)
    REAL(8) :: x, y, z, r, theta, phi
    INTEGER :: i
    CALL setup_dataset()
    CALL read_points(6, points)
    DO i = 1, SIZE(points, 2)
      CALL CAMBCORD(points(1,i), points(2,i), points(3,i), points(4,i), points(5,i), points(6,i), &
                    x, y, z, r, theta, phi, PESList(1))
      CALL put('C', (/ x, y, z, r, theta, phi /))
    END DO
  END SUBROUTINE

  SUBROUTINE run_spline1d()
    TYPE(Spline1D) :: s
    INTEGER :: n, ic1, icn, m, i
    REAL(8) :: ends(2), value, slope
    REAL(8), ALLOCATABLE :: t(:)
    READ(*, *) n, ic1, icn, m
    CALL read_hex(ends)
    CALL Allocate_Spline1D(s, n)
    CALL read_hex(s%X)
    CALL read_hex(s%F)
    ALLOCATE(t(m))
    CALL read_hex(t)
    CALL Compute_Spline1D(s, ic1, ends(1), icn, ends(2))
    CALL put('cm', s%Coefs)
    DO i = 1, m
      CALL SPLINT(s%X, s%F, s%Coefs, n, t(i), value, slope)
      CALL put('Q', (/ Evaluate_Spline1D(s, t(i)), Evaluate_Spline1D_Derivative(s, t(i)), value, slope /))
    END DO
  END SUBROUTINE

  SUBROUTINE run_spline2d()
    TYPE(Spline2D) :: s
    INTEGER :: nx, ny, m, i
    REAL(8) :: q(2), out(3)
    REAL(8), ALLOCATABLE :: f(:)
    LOGICAL :: outside
    READ(*, *) nx, ny, m
    CALL Allocate_Spline2D(s, nx, ny)
    CALL read_hex(s%X)
    CALL read_hex(s%Y)
    ALLOCATE(f(nx * ny))
    CALL read_hex(f)
    s%F = RESHAPE(f, (/ nx, ny /))
    CALL Compute_Spline2D(s)
    CALL put('cx', RESHAPE(s%CoefsX, (/ nx * ny /)))
    CALL put('cy', RESHAPE(s%CoefsY, (/ nx * ny /)))
    CALL put('cxy', RESHAPE(s%CoefsXY, (/ nx * ny /)))
    DO i = 1, m
      CALL read_hex(q)
      out = 0.0_8
      CALL Evaluate_Spline2D_Alt(s, q(1), q(2), out, outside)
      IF (outside) THEN
        PRINT '(A)', 'Q outside'
      ELSE
        CALL put('Q', out)
      END IF
    END DO
  END SUBROUTINE

  SUBROUTINE run_spline3d()
    TYPE(Spline3DPeriodic) :: s
    INTEGER :: nx, ny, nz, m, i
    REAL(8) :: q(3), f0, dx, dy, dz
    REAL(8), ALLOCATABLE :: f(:), d(:)
    READ(*, *) nx, ny, nz, m
    CALL Allocate_Spline3DPeriodic(s, nx, ny, nz)
    CALL read_hex(s%X)
    CALL read_hex(s%Y)
    CALL read_hex(s%Z)
    ALLOCATE(f(nx * ny * nz), d(nx * ny))
    CALL read_hex(f)
    s%F = RESHAPE(f, (/ nx, ny, nz /))
    CALL read_hex(d)
    s%Der1 = RESHAPE(d, (/ nx, ny /))
    CALL read_hex(d)
    s%Der2 = RESHAPE(d, (/ nx, ny /))
    CALL Compute_Spline3D_Periodic(s)
    CALL put('cx', RESHAPE(s%CoefsX, (/ nx * ny * nz /)))
    CALL put('cy', RESHAPE(s%CoefsY, (/ nx * ny * nz /)))
    CALL put('cxy', RESHAPE(s%CoefsXY, (/ nx * ny * nz /)))
    CALL put('cz', RESHAPE(s%CoefsZ, (/ nx * ny * nz /)))
    CALL put('czx', RESHAPE(s%CoefsZX, (/ nx * ny * nz /)))
    CALL put('czy', RESHAPE(s%CoefsZY, (/ nx * ny * nz /)))
    CALL put('czxy', RESHAPE(s%CoefsZXY, (/ nx * ny * nz /)))
    DO i = 1, m
      CALL read_hex(q)
      CALL Evaluate_Spline3D_Periodic(s, q(1), q(2), q(3), f0, dx, dy, dz)
      CALL put('Q', (/ f0, dx, dy, dz /))
    END DO
  END SUBROUTINE

  SUBROUTINE run_trispl()
    INTEGER :: n
    REAL(8), ALLOCATABLE :: x(:), y(:), cm(:)
    READ(*, *) n
    ALLOCATE(x(n), y(n), cm(n))
    CALL read_hex(x)
    CALL read_hex(y)
    cm = 0.0_8
    CALL TRISPL(n, x, y, cm)
    CALL put('cm', cm)
  END SUBROUTINE

  SUBROUTINE run_ubicar()
    INTEGER :: n, m, i, klo
    REAL(8), ALLOCATABLE :: x(:), t(:)
    LOGICAL :: outside
    READ(*, *) n, m
    ALLOCATE(x(n), t(m))
    CALL read_hex(x)
    CALL read_hex(t)
    DO i = 1, m
      klo = 0
      CALL UBICAR(x, t(i), n, klo, outside)
      WRITE(*, '(A,1X,I0,1X,I0,1X,L1)') 'U', i, klo, outside
    END DO
  END SUBROUTINE

END PROGRAM
