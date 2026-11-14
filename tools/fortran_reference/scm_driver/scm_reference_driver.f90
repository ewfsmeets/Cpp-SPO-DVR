! Reference driver for the C++ port of the SCM distorted-surface correction (src/potential/scm/).
!
! Calls the SCM and CRP routines of the Fortran reference build (objects of the reference install,
! see build.sh) in the order SPO-DVR uses for SystemID CRPDistrtSurface: Set_PES_DataSet and
! InitialiseDistPESModule (mode 1 = SCM; H2SurfaceSystem.f90 performActionH2SurfaceSystem), then
! GetCRPDistSurfParams (InputNew.f90 via setupParametersH2SurfaceSystem; it sets the CRP up lazily),
! InitializeCRPDistSurfModule (Main.f90 via newH2SurfaceSystem; the program's repeated
! Set_PES_DataSet and InitialiseDistPESModule there are no-ops and left out), then
! SetupDistPESModule (Potential.f90) once with RArray = the r of every point, so that point k uses
! the SCM parameters of "grid point" ir = k. The set-up can run only once per process
! (Setup_3Body_Rydberg_Par allocates without a guard): one case per run. Every run ends with the
! line SCM_DRIVER_DONE in the output file; bad arguments exit with status 2, failures with 1.
!
! Usage: scm_reference_driver <dataset> <params> <slab> <slab number> <cutoff> <T|F> <output> < points
!   dataset      the CRP data set directory (Set_PES_DataSet; keep it short, as for the CRP driver)
!   params slab  params.dat and the slab position file, at most 100 characters each
!                (InitialiseDistPESModule takes CHARACTER(100))
!   slab number  PosFileSlabNum; cutoff: SCMCutoff (bohr); T|F: EnableSurfExpan
! Points (stdin): the number of points, then one point per line, u v Z r theta phi as 16
! hexadecimal IEEE digits each (u, v: skewed lattice coordinates in bohr, as SPO-DVR passes them).
!
! Output (the file <output>): one line per double, "<tag> <index> <16 hex digits> <ES26.17E3>"
! (the CRP driver's form), other lines "<tag> <value>...". In this order:
!   rescaling, slab_number, cutoff                         the settings
!   point_group, lattice_x, lattice_y, z_asymptotic,       GetCRPDistSurfParams (lattice_x/_y are
!     skew_angle                                           the slab's a_T with rescaling)
!   ideal_lattice_constant, slab_lattice_constant, scale   a0 of the CRP, a_T of the slab file, and
!                                                          s = a_T/a0 with rescaling, else 1
!   cos_skew, sin_skew                                     the (u, v) -> (x, y) matrix elements
!   parameter_version, parameter_layer1 1..18,             params.dat as read (MySCMData)
!     parameter_layer2 1..18
!   slab_tag, slab_structure, slab_size (nx ny nz),        the slab file header as read
!     slab_configurations, slab_version, slab_max_atoms
!   atoms N, then per atom i: atom_layer i <layer>,        the atoms the module keeps (MySCMData;
!     distorted_x/_y/_z i, ideal_x/_y/_z i                 version 2: also those of layers > 2)
!   points N, then per point k: x k, y k                   the Cartesian centre of mass, as below
!     delta_v k                                            GetDistPESEnergy(x, y, Z, r, theta, phi)
!                                                          in Hartree: the correction alone
!     v_crp k                                              the CRP at (u/s, v/s, Z, r, theta, phi), eV
!     v_total k                                            GetCRPDistSurfPotential(u, v, ...), eV
!   SCM_DRIVER_DONE
! x and y are computed here as in DistortedSurfaceWrapper.f90 (MATMUL with the matrix of
! GetCRPDistSurfParams), u/s and v/s as in GetCRPDistSurfPotential; the driver checks that
! v_total = v_crp + delta_v / Constants_eV2Ha bitwise at every point, which holds only if both are
! the module's, and stops with status 1 otherwise. Module messages go to stdout.
PROGRAM scm_reference_driver
  USE DistortedSurfaceWrapper, ONLY: InitialiseDistPESModule, InitializeCRPDistSurfModule, &
                                     GetCRPDistSurfParams, SetupDistPESModule, GetDistPESEnergy, &
                                     GetCRPDistSurfPotential
  USE CRP_QD_Wrapper, ONLY: Set_PES_DataSet, Get_CRP_QD_Wrapper_Potential_6D
  USE SCM_Types, ONLY: MySCMData
  USE Constants, ONLY: Constants_eV2Ha
  USE ISO_FORTRAN_ENV, ONLY: OUTPUT_UNIT, ERROR_UNIT
  IMPLICIT NONE

  INTEGER, PARAMETER :: SCM_MODE = 1  ! DistPESMode of SPO-DVR: 1 = SCM
  CHARACTER(LEN=512) :: dataset, argument, output
  CHARACTER(LEN=100) :: params_file, slab_file
  CHARACTER(LEN=10)  :: point_group
  LOGICAL :: rescaling
  INTEGER :: slab_number, out, k, i, mismatches, status
  REAL(8) :: cutoff, xlat, ylat, zasymptotic, skew, scale
  REAL(8) :: to_cartesian(2, 2), xy(2), uv_ideal(2), delta_v, v_crp, v_total
  REAL(8), ALLOCATABLE :: points(:, :)

  IF (COMMAND_ARGUMENT_COUNT() /= 7) CALL usage('7 arguments expected')
  CALL GET_COMMAND_ARGUMENT(1, dataset)
  params_file = short_path(2)
  slab_file = short_path(3)
  CALL GET_COMMAND_ARGUMENT(4, argument)
  READ(argument, *, IOSTAT=status) slab_number
  IF (status /= 0) CALL usage('the slab number must be an integer')
  CALL GET_COMMAND_ARGUMENT(5, argument)
  READ(argument, *, IOSTAT=status) cutoff
  IF (status /= 0) CALL usage('the cutoff must be a number')
  CALL GET_COMMAND_ARGUMENT(6, argument)
  rescaling = (TRIM(argument) == 'T')
  IF (.NOT. rescaling .AND. TRIM(argument) /= 'F') CALL usage('the rescaling argument must be T or F')
  CALL GET_COMMAND_ARGUMENT(7, output)
  CALL read_points(points)

  ! The set-up, in SPO-DVR's order.
  CALL Set_PES_DataSet(TRIM(dataset))
  IF (.NOT. InitialiseDistPESModule(params_file, slab_file, rescaling, cutoff, slab_number, SCM_MODE)) &
    CALL fail('InitialiseDistPESModule failed (see the messages above); check the slab file and the cutoff')
  CALL GetCRPDistSurfParams(xlat, ylat, zasymptotic, skew, point_group)
  CALL InitializeCRPDistSurfModule()
  IF (.NOT. SetupDistPESModule(points(4, :), SCM_MODE)) &
    CALL fail('SetupDistPESModule failed (see the messages above); check the params and slab files and the slab number')

  ! The (u, v) -> (x, y) matrix and the rescaling factor, as DistortedSurfaceWrapper.f90 has them.
  to_cartesian(1, :) = (/ 1.0, COS(skew) /)
  to_cartesian(2, :) = (/ 0.0, SIN(skew) /)
  scale = 1.
  IF (rescaling) scale = MySCMData%LatExpansion

  OPEN(NEWUNIT=out, FILE=TRIM(output), STATUS='REPLACE', ACTION='WRITE')
  WRITE(out, '(A,1X,L1)') 'rescaling', rescaling
  WRITE(out, '(A,1X,I0)') 'slab_number', slab_number
  CALL put('cutoff', 1, cutoff)
  WRITE(out, '(A,1X,A)') 'point_group', TRIM(point_group)
  CALL put('lattice_x', 1, xlat)
  CALL put('lattice_y', 1, ylat)
  CALL put('z_asymptotic', 1, zasymptotic)
  CALL put('skew_angle', 1, skew)
  CALL put('ideal_lattice_constant', 1, MySCMData%IdLatticeConstant)
  CALL put('slab_lattice_constant', 1, MySCMData%DistLatticeConstant)
  CALL put('scale', 1, scale)
  CALL put('cos_skew', 1, to_cartesian(1, 2))
  CALL put('sin_skew', 1, to_cartesian(2, 2))

  WRITE(out, '(A,1X,I0)') 'parameter_version', MySCMData%ParameterVersion
  DO k = 1, 18
    CALL put('parameter_layer1', k, MySCMData%RydbergParams(1)%P(k))
  END DO
  DO k = 1, 18
    CALL put('parameter_layer2', k, MySCMData%RydbergParams(2)%P(k))
  END DO

  WRITE(out, '(A,1X,A)') 'slab_tag', TRIM(MySCMData%SlabFileTag)
  WRITE(out, '(A,1X,A)') 'slab_structure', TRIM(MySCMData%SlabSurfStruct)
  WRITE(out, '(A,3(1X,I0))') 'slab_size', MySCMData%SlabAtomSize
  WRITE(out, '(A,1X,I0)') 'slab_configurations', MySCMData%NumSlabConfigs
  WRITE(out, '(A,1X,I0)') 'slab_version', MySCMData%PosFileVersion
  WRITE(out, '(A,1X,I0)') 'slab_max_atoms', MySCMData%MaxNumSurfAtoms
  WRITE(out, '(A,1X,I0)') 'atoms', MySCMData%NumVectors
  DO i = 1, MySCMData%NumVectors
    WRITE(out, '(A,1X,I0,1X,I0)') 'atom_layer', i, MySCMData%AtomType(i)
    CALL put('distorted_x', i, MySCMData%VectorsDisplaced(3 * i - 2))
    CALL put('distorted_y', i, MySCMData%VectorsDisplaced(3 * i - 1))
    CALL put('distorted_z', i, MySCMData%VectorsDisplaced(3 * i))
    CALL put('ideal_x', i, MySCMData%VectorsIdeal(3 * i - 2))
    CALL put('ideal_y', i, MySCMData%VectorsIdeal(3 * i - 1))
    CALL put('ideal_z', i, MySCMData%VectorsIdeal(3 * i))
  END DO

  ! The points: the correction alone at the Cartesian centre of mass, the CRP at the ideal-lattice
  ! (u, v), and SPO-DVR's total (the CRP index arguments iZ, iR are unused there: 0, k).
  WRITE(out, '(A,1X,I0)') 'points', SIZE(points, 2)
  mismatches = 0
  DO k = 1, SIZE(points, 2)
    xy = MATMUL(to_cartesian, points(1:2, k))
    delta_v = GetDistPESEnergy(xy(1), xy(2), points(3, k), points(4, k), points(5, k), points(6, k), k, SCM_MODE)
    uv_ideal = points(1:2, k) / scale
    v_crp = Get_CRP_QD_Wrapper_Potential_6D(uv_ideal(1), uv_ideal(2), points(3, k), points(4, k), &
                                            points(5, k), points(6, k), 0, k)
    v_total = GetCRPDistSurfPotential(points(1, k), points(2, k), points(3, k), points(4, k), &
                                      points(5, k), points(6, k), 0, k)
    IF (TRANSFER(v_total, 0_8) /= TRANSFER(v_crp + delta_v / Constants_eV2Ha, 0_8)) mismatches = mismatches + 1
    CALL put('x', k, xy(1))
    CALL put('y', k, xy(2))
    CALL put('delta_v', k, delta_v)
    CALL put('v_crp', k, v_crp)
    CALL put('v_total', k, v_total)
  END DO
  IF (mismatches > 0) THEN
    WRITE(out, '(A,1X,I0,A)') 'error:', mismatches, ' points with v_total /= v_crp + delta_v / Constants_eV2Ha'
    CLOSE(out)
    CALL fail('v_total differs from v_crp + delta_v / Constants_eV2Ha: the driver''s x, y or u/s, v/s ' // &
              'are not the module''s; make them as DistortedSurfaceWrapper.f90 does')
  END IF
  WRITE(out, '(A)') 'SCM_DRIVER_DONE'
  CLOSE(out)

CONTAINS

  ! Stops with status 1 after printing the message on stderr (a plain STOP would exit with 0),
  ! after the module messages already printed on stdout.
  SUBROUTINE fail(message)
    CHARACTER(LEN=*), INTENT(IN) :: message
    FLUSH(OUTPUT_UNIT)
    WRITE(ERROR_UNIT, '(A)') 'error: ' // message
    FLUSH(ERROR_UNIT)
    STOP 1
  END SUBROUTINE

  ! A bad argument: what is wrong and the usage line on stderr, exit status 2.
  SUBROUTINE usage(problem)
    CHARACTER(LEN=*), INTENT(IN) :: problem
    WRITE(ERROR_UNIT, '(A)') 'scm_reference_driver: ' // problem
    WRITE(ERROR_UNIT, '(A)') 'usage: scm_reference_driver <dataset> <params> <slab> <slab number> <cutoff> ' // &
                             '<T|F> <output> < points'
    FLUSH(ERROR_UNIT)
    STOP 2
  END SUBROUTINE

  ! Command-line argument n as a CHARACTER(100) path; a longer one, which would be cut, is a bad
  ! argument.
  FUNCTION short_path(n) RESULT(path)
    INTEGER, INTENT(IN) :: n
    CHARACTER(LEN=100) :: path
    CHARACTER(LEN=512) :: text
    CALL GET_COMMAND_ARGUMENT(n, text)
    IF (LEN_TRIM(text) > LEN(path)) CALL usage('a file path longer than 100 characters: ' // TRIM(text))
    path = text(1:LEN(path))
  END FUNCTION

  ! A double from 16 hex digits (read as two 32-bit halves: no overflow of INTEGER(8)).
  REAL(8) FUNCTION from_hex(token)
    CHARACTER(LEN=*), INTENT(IN) :: token
    INTEGER(8) :: high, low
    READ(token(1:8), '(Z8)') high
    READ(token(9:16), '(Z8)') low
    from_hex = TRANSFER(IOR(ISHFT(high, 32), low), 1.0_8)
  END FUNCTION

  ! The point list from stdin: points(:, k) = u v Z r theta phi of point k.
  SUBROUTINE read_points(values)
    REAL(8), ALLOCATABLE, INTENT(OUT) :: values(:, :)
    CHARACTER(LEN=16) :: tokens(6)
    INTEGER :: n, p, c
    READ(*, *) n
    ALLOCATE(values(6, n))
    DO p = 1, n
      READ(*, *) tokens
      DO c = 1, 6
        values(c, p) = from_hex(tokens(c))
      END DO
    END DO
  END SUBROUTINE

  SUBROUTINE put(tag, index, value)
    CHARACTER(LEN=*), INTENT(IN) :: tag
    INTEGER, INTENT(IN) :: index
    REAL(8), INTENT(IN) :: value
    WRITE(out, '(A,1X,I0,1X,Z16.16,1X,ES26.17E3)') tag, index, TRANSFER(value, 0_8), value
  END SUBROUTINE

END PROGRAM scm_reference_driver
