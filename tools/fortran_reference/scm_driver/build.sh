#!/usr/bin/env bash
# Builds the SCM reference driver (scm_reference_driver.f90) in a scratch directory, linked against
# the SCM and CRP library objects of the reference install (reference_config.sh: OBJS/, checked
# against its SHA256SUMS), so that the driver runs exactly the code generated for the reference
# executable. The .mod files there are specific to the gfortran version that built them, so the
# installed gfortran must be the recorded one (BUILDINFO.tsv package.build.gfortran-15).
#
# Usage: tools/fortran_reference/scm_driver/build.sh [build_dir]
#   build_dir  a directory that does not exist yet (default: a new mktemp -d)
# Prints the path of the executable on the last line.
# Environment: SPODVR_FORTRAN_REF_DIR, SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT (reference_config.sh).
# Run time: seconds. Docs: tools/fortran_reference/README.md ("The SCM driver").
set -euo pipefail
source_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${source_dir}/../reference_config.sh"

(($# <= 1)) || { echo "usage: scm_driver/build.sh [build_dir]" >&2; exit 2; }
# The objects the call chain needs (each one is: the link fails without it): the SCM modules
# (DistortedSurfaceWrapper -> SCM_Wrapper -> SCM_IO, SCM_Potential, SCM_Types) and the CRP library
# as in crp_driver/build.sh (errortrap: the 5.24 accessor Get_CRP_QD_Wrapper_Additional_Parameters,
# which GetCRPDistSurfParams calls, uses AbortWithError). constants.o and CRP_Constants.o are left
# out: no linked object references a symbol of theirs (nm -u), only PARAMETERs of their .mod files.
build_fortran_driver "${1:-}" scm_reference_driver "${source_dir}/scm_reference_driver.f90" -- \
  DistortedSurfaceWrapper SCM_Wrapper SCM_IO SCM_Potential SCM_Types \
  CRP_Wrapper CRP_QD_Wrapper CRP_General CRP_GasPhase \
  CRP_C3v_AtomSurface CRP_C3v_MoleculeSurface CRP_C4v_AtomSurface CRP_C4v_MoleculeSurface \
  CRP_C6v_AtomSurface CRP_C6v_MoleculeSurface Spline UnitConversion MPI_wrapper errortrap
