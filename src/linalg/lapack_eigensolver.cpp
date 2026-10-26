// LAPACK (dsyevd) backend of the symmetric eigensolver; C-API boundary to LAPACK.
#include <cstddef>
#include <stdexcept>
#include <string>

#include "linalg/symmetric_eigensolver.hpp"

// This file is the C-API boundary to LAPACK: raw pointers are used only here.
extern "C" void dsyevd_(const char* jobz, const char* uplo, const int* n, double* a, const int* lda,
                        double* w, double* work, const int* lwork, int* iwork, const int* liwork,
                        int* info, std::size_t jobz_len, std::size_t uplo_len);

namespace spodvr {

EigenDecomposition LapackEigensolver::solve(std::span<const double> matrix, int n) const {
  if (n < 1 || matrix.size() < static_cast<std::size_t>(n) * static_cast<std::size_t>(n))
    throw std::invalid_argument("LapackEigensolver: invalid matrix size");

  EigenDecomposition result;
  result.n = n;
  result.values.resize(static_cast<std::size_t>(n));
  result.vectors.assign(matrix.begin(), matrix.begin() + static_cast<std::ptrdiff_t>(n) * n);

  const char jobz = 'V';
  const char uplo = 'L';
  int info = 0;

  // Workspace query.
  int lwork = -1;
  int liwork = -1;
  double work_size = 0.0;
  int iwork_size = 0;
  dsyevd_(&jobz, &uplo, &n, result.vectors.data(), &n, result.values.data(), &work_size, &lwork,
          &iwork_size, &liwork, &info, 1, 1);
  if (info != 0) throw std::runtime_error("LapackEigensolver: dsyevd workspace query failed");

  lwork = static_cast<int>(work_size);
  liwork = iwork_size;
  std::vector<double> work(static_cast<std::size_t>(lwork));
  std::vector<int> iwork(static_cast<std::size_t>(liwork));
  dsyevd_(&jobz, &uplo, &n, result.vectors.data(), &n, result.values.data(), work.data(), &lwork,
          iwork.data(), &liwork, &info, 1, 1);
  if (info != 0)
    throw std::runtime_error("LapackEigensolver: dsyevd failed (info = " + std::to_string(info) + ")");
  return result;
}

}  // namespace spodvr
