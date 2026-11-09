/* Prints the Fourier coefficient tables of the Fortran CRP library bit for bit.
 * The tables are PRIVATE module variables (filled by Initialize_CRP_Wrapper_Module), but gfortran
 * emits them as global symbols, so they can be read from C. Output lines:
 *   <table> <index (Fortran column-major order, 0-based)> <hex bits> <value>
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

extern double __crp_c6v_moleculesurface_MOD_fouri6d[144];
extern double __crp_c3v_moleculesurface_MOD_fouri6d[144];
extern double __crp_c6v_atomsurface_MOD_fouri3d[49];
extern double __crp_c3v_atomsurface_MOD_fouri3d[100];

static void dump(const char* name, const double* table, int n) {
  for (int i = 0; i < n; ++i) {
    uint64_t bits;
    memcpy(&bits, &table[i], sizeof bits);
    printf("%s %d %016llx %.17e\n", name, i, (unsigned long long)bits, table[i]);
  }
}

void dump_crp_tables(void) {
  dump("fouri6d_c6v", __crp_c6v_moleculesurface_MOD_fouri6d, 144);
  dump("fouri6d_c3v", __crp_c3v_moleculesurface_MOD_fouri6d, 144);
  dump("fouri3d_c6v", __crp_c6v_atomsurface_MOD_fouri3d, 49);
  dump("fouri3d_c3v", __crp_c3v_atomsurface_MOD_fouri3d, 100);
  fflush(stdout);
}
