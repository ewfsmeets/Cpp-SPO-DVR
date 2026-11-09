// Fourier coefficient tables of the CRP interpolation. Fortran: FOUR6D (CRP_C6v/C3v_MoleculeSurface),
// FOUR3D (CRP_C6v_AtomSurface, CRP_C3v_AtomSurface). The DATA lists are written out with their
// repeat counts expanded, column by column (Fortran storage order).
#include "crplib/fourier_tables.hpp"

#include "crplib/constants.hpp"

namespace crplib::detail {

namespace {

// Fortran: FOUR6D DATA tour (zer, un, mun, deu, mdeu = 0, 1, -1, 2, -2)
constexpr std::array<double, 144> tour_6d = {
    1, 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   // 12*un
    2, -2, -1, 1,  -2, 2,  -1, 1,  -1, 1,  -1, 1,   // deu, mdeu, mun, un, mdeu, deu, mun, un, mun, un, mun, un
    2, 2,  -1, 1,  -2, -2, 1,  -1, -1, 1,  1,  -1,  // 2*deu, mun, un, 2*mdeu, un, 2*mun, 2*un, mun
    2, -2, -1, -1, 2,  -2, 1,  1,  -1, -1, 1,  1,   // deu, mdeu, 2*mun, deu, mdeu, 2*un, 2*mun, 2*un
    0, 0,  1,  1,  0,  0,  -1, 1,  -1, -1, 1,  -1,  // 2*zer, 2*un, 2*zer, mun, un, 2*mun, un, mun
    0, 0,  1,  1,  0,  0,  1,  -1, -1, -1, -1, 1,   // 2*zer, 2*un, 2*zer, un, 4*mun, un
    0, 0,  -1, 1,  0,  0,  -1, -1, 1,  -1, 1,  1,   // 2*zer, mun, un, 2*zer, 2*mun, un, mun, 2*un
    1, 1,  1,  -1, -1, -1, -1, 1,  1,  -1, -1, 1,   // 3*un, 4*mun, 2*un, 2*mun, un
    1, -1, 1,  -1, -1, 1,  1,  -1, 1,  -1, 1,  -1,  // un, mun, un, 2*mun, 2*un, mun, un, mun, un, mun
    1, -1, 1,  1,  1,  -1, -1, -1, 1,  1,  -1, -1,  // un, mun, 3*un, 3*mun, 2*un, 2*mun
    2, 2,  -1, -1, 2,  2,  -1, -1, -1, -1, -1, -1,  // 2*deu, 2*mun, 2*deu, 6*mun
    0, 0,  1,  -1, 0,  0,  -1, -1, -1, 1,  1,  1};  // 2*zer, un, mun, 2*zer, 3*mun, 3*un

// Fortran: FOUR6D, fouri6D(i, j) = 3.0 * r2 * tour (j = 1..4, 8..10), r1 * tour (5..7),
// r2 * tour (11), r1 / 3.0 * tour (12), with r1 = usq3 / 4.00 and r2 = 1.00 / 36.00.
constexpr std::array<double, 144> make_table_6d() {
  constexpr double r1 = usq3 / 4.00;
  constexpr double r2 = 1.00 / 36.00;
  std::array<double, 144> table{};
  for (int j = 0; j < 12; ++j) {
    double factor = 3.0 * r2;
    if (j >= 4 && j <= 6) factor = r1;
    if (j == 10) factor = r2;
    if (j == 11) factor = r1 / 3.0;
    for (int i = 0; i < 12; ++i) table[i + 12 * j] = factor * tour_6d[i + 12 * j];
  }
  return table;
}

// Fortran: CRP_C6v_AtomSurface.f90 FOUR3D DATA tour (7 x 7)
constexpr std::array<double, 49> tour_3d_c6v = {1.0,
                                                3.0,
                                                2.0,
                                                6.0,
                                                12.0,
                                                6.0,
                                                6.0,   //
                                                2.0,
                                                -2.0,
                                                -2.0,
                                                2.0,
                                                -8.0,
                                                8.0,
                                                0.0,   //
                                                2.0,
                                                -2.0,
                                                4.0,
                                                -4.0,
                                                4.0,
                                                2.0,
                                                -6.0,  //
                                                2.0,
                                                6.0,
                                                -2.0,
                                                -6.0,
                                                0.0,
                                                0.0,
                                                0.0,   // 2.0,6.0,-2.0,-6.0,3*0.0
                                                2.0,
                                                -2.0,
                                                -2.0,
                                                2.0,
                                                4.0,
                                                -4.0,
                                                0.0,   //
                                                1.0,
                                                -1.0,
                                                2.0,
                                                -2.0,
                                                -4.0,
                                                -2.0,
                                                6.0,   //
                                                0.666666666666667,
                                                2.0,
                                                1.333333333333,
                                                4.0,
                                                -4.0,
                                                -2.0,
                                                -2.0};

// Fortran: CRP_C3v_AtomSurface.f90 FOUR3D DATA tour (10 x 10)
constexpr double a = 5.1961524227066318805600;
constexpr double b = 20.7846096908265275222400;
constexpr double c = 15.5884572681198956416800;
constexpr double d = 10.3923048454132637611200;
constexpr std::array<double, 100> tour_3d_c3v = {
    3.00, 9.00,  3.00,  9.00,  18.00,  18.00,  18.00,  9.00,  18.00,  3.00,   // ..., 3*18.00, ...
    6.00, -6.00, -3.00, 3.00,  -12.00, 24.00,  0.00,   3.00,  -12.00, -3.00,  //
    0.00, 0.00,  -a,    -a,    -b,     0.00,   0.00,   a,     b,      a,      // 2*0.00, 2*-a, -b, 2*0.00, a, b, a
    6.00, -6.00, 6.00,  -6.00, 6.00,   6.00,   -18.00, -6.00, 6.00,   6.00,   //
    6.00, 18.00, -3.00, -9.00, 0.00,   0.00,   0.00,   -9.00, 0.00,   -3.00,  //
    0.00, 0.00,  a,     -c,    0.00,   0.00,   0.00,   c,     0.00,   -a,     // 2*0.00, a, -c, 3*0.00, c, 0.00, -a
    6.00, -6.00, -3.00, 3.00,  6.00,   -12.00, 0.00,   3.00,  6.00,   -3.00,  //
    0.00, 0.00,  -a,    -a,    d,      0.00,   0.00,   a,     -d,     a,      // 2*0.00, 2*-a, d, 2*0.00, a, -d, a
    3.00, -3.00, 3.00,  -3.00, -6.00,  -6.00,  18.00,  -3.00, -6.00,  3.00,   //
    2.00, 6.00,  2.00,  6.00,  -6.00,  -6.00,  -6.00,  6.00,  -6.00,  2.00};  // ..., 3*-6.00, ...

template <std::size_t N>
constexpr std::array<double, N> divided(const std::array<double, N>& tour, double divisor) {
  std::array<double, N> table{};
  for (std::size_t i = 0; i < N; ++i) table[i] = tour[i] / divisor;
  return table;
}

constexpr std::array<double, 144> table_6d = make_table_6d();
constexpr std::array<double, 49> table_3d_c6v = divided(tour_3d_c6v, 36.00);    // Fortran: tour / 36.00
constexpr std::array<double, 100> table_3d_c3v = divided(tour_3d_c3v, 108.00);  // Fortran: tour / 108.00

}  // namespace

const std::array<double, 144>& fourier_table_molecule_hexagonal() { return table_6d; }
const std::array<double, 49>& fourier_table_atom_c6v() { return table_3d_c6v; }
const std::array<double, 100>& fourier_table_atom_c3v() { return table_3d_c3v; }

}  // namespace crplib::detail
