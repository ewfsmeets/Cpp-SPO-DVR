! Validation driver (angular matrices, docs/validation_thesis_fig5_3_dz): the Fortran angular grid
! and DVR<->FBR transform of the thesis grid (JMax = MjMax = 28, homonuclear, JIn = MjIn = 0, nphi a
! multiple of 6), set up as setupdvr.f90 SetupDVRProj does it. Linked against the object files of
! the reference executable (~/work/SPO-DVR_v5.19/procfiles.dir/OBJS.dir, read only); see angular_matrices.py.
! Writes into the current directory, in the layout of angular_check.cpp:
!   thetas.bin, weights.bin, matrix.bin (FBRToDVRMatrices), basis.txt, roundtrip.bin
PROGRAM AngularFortran
USE Constants
USE Legendre
USE FFTs
USE AngularGrid
USE Transform
IMPLICIT NONE
INTEGER, PARAMETER :: JMax = 28, MjMax = 28, JIn = 0, MjIn = 0, PhiMultiple = 6
TYPE(AngularGridType)         :: Grid
TYPE(GenComplexTransformType) :: Tr
INTEGER                       :: nPhi, n, i, u
REAL                          :: dPhi
COMPLEX, ALLOCATABLE          :: Fbr(:,:), Dvr(:,:,:)

    ! as SetupDVRProj (SymmetryInPhi = .FALSE.)
    CALL SetupConstantsModule()
    CALL SetupFFTModule
    CALL SetupLegendreModule( JMax + 1 )
    CALL SetupAngularGridModule( JMax )
    nPhi = GetNearestFFTSize( 2 * MjMax + 1, PhiMultiple )
    dPhi = Constants_2PI / nPhi
    CALL InitAngularGrid( Grid, .TRUE., .FALSE., nPhi, dPhi, JMax, JIn, MjIn )
    CALL SetupTransformModule( JMax )
    CALL InitGenComplexTransform( Grid, Tr, MjMax, PhiMultiple )
    n = Tr%NrOfRotIndices

    OPEN( NEWUNIT = u, FILE = 'thetas.bin', ACCESS = 'STREAM', FORM = 'UNFORMATTED', STATUS = 'REPLACE' )
    WRITE( u ) Grid%TheThetas
    CLOSE( u )
    OPEN( NEWUNIT = u, FILE = 'weights.bin', ACCESS = 'STREAM', FORM = 'UNFORMATTED', STATUS = 'REPLACE' )
    WRITE( u ) Grid%TheWeights
    CLOSE( u )
    OPEN( NEWUNIT = u, FILE = 'matrix.bin', ACCESS = 'STREAM', FORM = 'UNFORMATTED', STATUS = 'REPLACE' )
    WRITE( u ) Tr%FBRToDVRMatrices
    CLOSE( u )
    OPEN( NEWUNIT = u, FILE = 'basis.txt', STATUS = 'REPLACE' )
    WRITE( u, '(I0,1X,I0,1X,I0)' ) n, Grid%nPhis, Grid%nThetas
    DO i = 1, n
       WRITE( u, '(I0,1X,I0)' ) Tr%JBasisSet( i ), Tr%MjBasisSet( i )
    END DO
    CLOSE( u )

    ! the round trip DVR->FBR(FBR->DVR(e_J)) of every mJ = 0 basis function (one column: the
    ! Fortran transforms every column separately)
    ALLOCATE( Fbr( n, 1 ), Dvr( nPhi + 1, Grid%nThetas, 1 ) )
    OPEN( NEWUNIT = u, FILE = 'roundtrip.bin', ACCESS = 'STREAM', FORM = 'UNFORMATTED', STATUS = 'REPLACE' )
    DO i = 1, n
       IF ( Tr%MjBasisSet( i ) /= 0 ) CYCLE
       Fbr = CMPLX( 0.0, 0.0 )
       Fbr( i, 1 ) = CMPLX( 1.0, 0.0 )
       Dvr = CMPLX( 0.0, 0.0 )
       CALL FBRToDVRGenComplexTransform( Tr, Fbr, Dvr, 1 )
       CALL DVRToFBRGenComplexTransform( Tr, Dvr, Fbr, 1 )
       WRITE( u ) Fbr( :, 1 )
    END DO
    CLOSE( u )
    PRINT '(A,I0,A,I0,A,I0)', 'angular_fortran: ', n, ' basis functions, ntheta ', Grid%nThetas, ', nphi ', nPhi
END PROGRAM AngularFortran
