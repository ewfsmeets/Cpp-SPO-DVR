#pragma once
// Errors of the CRP library. The Fortran library ends the program with STOP (exit status 0) where
// these are thrown.

#include <stdexcept>

namespace crplib {

// The data set cannot be used: a file is missing, malformed or inconsistent. The message names
// the file (and the line, for reading errors).
class DataError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

// The PES cannot be evaluated at this geometry (Fortran: STOP), e.g. an atom outside the z range
// of the atom-surface data. The message names the geometry.
class EvaluationError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

}  // namespace crplib
