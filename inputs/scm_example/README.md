# SCM technical example

**What it is.** A small, complete input for the SCM distorted-surface potential (SystemID
`CRPDistrtSurface`, [docs/scm.md](../../docs/scm.md)): the SRP48 H₂/Cu(111) CRP data set with the
example parameter and position files of `potentials/H2Cu111SRP_C3v.dir/` (lattice rescaling on),
H₂ (v, j) = (0, 0) at 0.5–0.9 eV, 3000 steps to t = 7500 a.u. It is the calculation of the gate case
`tests/reference/scm_cu111_srp48_small` (only the title differs), which compares it with the Fortran
in gates 1–4 ([docs/validation.md](../../docs/validation.md)); that directory's `DVR.inp` is the
Fortran input. Use it to see the SCM keys and the set-up report, and as a template for method tests.

**What it is not: physics.**

- One snapshot, the only one in the example position file. An SCM result is an average over many
  snapshots ([ensembles](../../docs/scm.md#ensembles-by-hand)).
- Test-sized grids (nX = nY = 6, JMax 4, 3 θ points): the run warns that the rotational and the
  diffraction basis are too small for all open channels, and the x, y grid cannot converge across
  the jump of V at the cell edges (SCM D10).
- The model itself has open questions (the units of the bond-length clamp, the cell-edge jump, the
  holes below the CRP's Zmin, …): [docs/scm_open_questions.md](../../docs/scm_open_questions.md).

Do not read its probabilities as an SCM result.

**Running.** From the repository root (the data paths in `input.json` are relative to it):

```bash
./build/spodvr inputs/scm_example/input.json --threads 4 --output-dir runs/scm_example
```

This took 1 min 34 s on 4 threads of an AMD Ryzen 7 5800X (2026-11-15). For the set-up report
alone, add `"Potential": true` to a copy of the input (the 6D potential took about 4 s). The block
"SCM DISTORTED SURFACE" of the output is explained in
[docs/scm.md](../../docs/scm.md#the-scm-set-up-report); the run's `analysis.json` holds its SCM
identity (`scm`, [docs/output_files.md](../../docs/output_files.md)).
