# example4/ — SEAS Benchmark BP8-QD-GS

Gaussian-source variant of BP8: quasi-dynamic rate-and-state slip on a
planar fault in an elastic whole-space, driven by pore-pressure diffusion
from a point injection at the fault center. See
`docs/SEAS_BP8_Benchmark_Description.pdf` for the full specification.

## Directory layout

```
example4/
  make_bp8_mesh.py     # generates triangular_mesh.gts (flat, uniform grid)
  parameter1.txt        # BP8 Table 1 parameters, see below for the format
  stiffness/             # triangular_mesh.gts + calc_nikkhoo_fs output
    triangular_mesh.gts
    trigreen_{22,23,32,33}_<rank>.bin
    position.bin
  out/                   # bp8_main's output (created empty; see .gitkeep)
```

## Mesh: exactly Omega_f, no locked buffer

Unlike example2/example3 (BP5/BP6, half-space, need a locked buffer around
the seismogenic patch to represent the surrounding fault), this mesh is
**exactly** the frictional domain Omega_f = (-lf,lf) x (-lf,lf) = 800m x
800m at the required 10m cell size (80x80 cells, 12,800 triangles). No
buffer is needed: BP8 Eq. (13) forces V=0 identically for all time outside
Omega_f, and a region with zero slip for all time contributes zero elastic
stress to its neighbors regardless of how far a mesh extends (a standard
property of boundary-integral crack/fault problems) — so meshing a locked
buffer around Omega_f would give the identical answer at a much higher
cost. This is *not* the same situation as BP5/6's transition zones, where
the "locked" region actually creeps at a nonzero (if slow) rate and does
contribute.

That said, this is an assumption worth checking empirically, and the
benchmark's own tip (Section 6) is to verify results are independent of
computational domain size — `make_bp8_mesh.py --domain <bigger>` makes
that straightforward if you want to confirm it (see "Domain-independence
check" below).

## Build

```bash
cd NikkhooWalter2015 && make -f Makefile.fs        # builds calc_nikkhoo_fs
cd ../src && make -f Makefile.bp8                   # builds bp8_main
```

## Run

```bash
# 1. Generate the mesh (once; already done for the default 800m/10m case)
cd example4/stiffness
python3 ../make_bp8_mesh.py --domain 800 --cell-size 10 --out triangular_mesh.gts

# 2. Compute the stiffness matrices (must use the SAME -np as step 3)
mpirun -np 4 ../../NikkhooWalter2015/calc_nikkhoo_fs

# 3. Run the simulation (same -np as step 2 -- both use an identical
#    element decomposition formula, so trigreen_*_<rank>.bin chunks must
#    line up with this run's per-rank element ranges)
cd ..
mpirun -np 4 ../src/bp8_main
```

`-np 4` with `OMP_NUM_THREADS=2` (8 total threads) is a reasonable default
for an 8-core machine; adjust both consistently for your hardware.

**Measured timing (this repo's 8-core dev machine, -np 4, OMP_NUM_THREADS=2):**
- Stiffness calculation (step 2, 12,800^2 element pairs x 2 full-space
  Green's function evaluations each): **~62 minutes**. The output
  (`trigreen_*.bin`, `position.bin`, ~4.9GB total for the 800m/10m mesh)
  is **not** committed to git (see `.gitignore` — this is regenerable data,
  and 4.9GB is far too large for a git repo) — you need to run step 2
  yourself before step 3. `triangular_mesh.gts` (500KB) *is* committed,
  so you don't need to regenerate the mesh itself unless you want a
  different resolution/domain size.
- Mechanical time-stepping run (step 3): the pore-pressure precompute
  (721 hourly snapshots on an 81x81 grid) takes well under a minute. The
  adaptive RK integration is the expensive part: **observed ~400
  simulated-seconds per wall-clock minute** early in the run (first ~19
  of 720 simulated hours, still well within the 100-hour injection
  window) — extrapolating, **the full 30-day run needs on the order of
  4-5 wall-days** on hardware like this. This will likely change (faster
  or slower) once slip actually accelerates past the early near-steady
  regime, so treat this as a rough planning number, not a guarantee.
  Budget accordingly, or run on more cores / a cluster.
- Every output file is flushed after each write (`write_all_output` in
  `3dtri_BP8.f90`), so an interrupted job (killed, preempted, crashed)
  keeps whatever it had written up to the last completed output step —
  there's no separate checkpoint/restart mechanism, so a genuinely
  interrupted run has to restart from t=0, but you won't lose *visibility*
  into how far it got.

## parameter1.txt format

Plain-text, one field group per line (this is a new, BP8-specific format —
much simpler than BP5/6's, since BP8's friction parameters are spatially
uniform, Table 1, rather than read from a var-*.dat profile file):

| Line | Contents |
|---|---|
| 1 | `foldername` — output directory (must exist, trailing `/`) |
| 2 | `stiffname` — path to calc_nikkhoo_fs output (trailing `/`) |
| 3 | `n_side dz_cell` — mesh cells across the full domain, cell size (m) |
| 4 | `mu nu rho` — shear modulus (Pa), Poisson's ratio, density (kg/m^3) |
| 5 | `seff0 tauinit` — initial effective normal stress, initial shear stress (Pa) |
| 6 | `a b Drs Vstar fstar` — rate-and-state friction parameters (Drs, Vstar in m, m/s) |
| 7 | `Vinit Vzero` — initial slip rate, floor value to avoid log(0) (m/s) |
| 8 | `Q0 Lfwid toff Lgauss` — injection rate (m^3/s), fault thickness (m), injection duration (s), Gaussian source width (m) |
| 9 | `alpha beta phi` — hydraulic diffusivity (m^2/s), compressibility (Pa^-1), porosity |
| 10 | `tf_end tint_out` — total simulation time, output interval (s) |

All values in SI base units (Pa, m, s, m/s) — the values in
`example4/parameter1.txt` are BP8 Table 1 converted to these units.

## Output format

Matches BP8 Section 4: 9 station time series (`fltst_strk±NNNdp±NNN`),
`global.dat` (Vmax, moment rate), and 10 spatial profile files
(`{slip,shear_stress}_{2,3}_{strike,depth}.dat`, `pore_pressure_{strike,depth}.dat`).

**One judgment call worth flagging**: the benchmark PDF's extracted text
renders compound field names with spaces ("slip 2", "shear stress 2") in
both the field-description table and the literal field-list line of its
own example files. Every other SEAS benchmark (BP1-BP7) uses
underscore_separated names in this exact format ("slip_2",
"shear_stress_2"), and the PDF's field table lists these as single
identifiers rather than prose, so `src/3dtri_BP8.f90` writes underscored
names, on the assumption the PDF's text extraction dropped underscores (a
known PDF-extraction artifact, especially in monospaced/code-formatted
text). If the actual CRESCENT DET upload server rejects the field-list
line, it's a single string constant to change (`FIELDS_TS` in
`open_output_files`, `src/3dtri_BP8.f90`).

## Domain-independence check

To verify the "mesh = exactly Omega_f" assumption (see above) holds in
practice, not just in principle, rerun steps 1-3 with a larger domain,
e.g.:

```bash
python3 make_bp8_mesh.py --domain 1200 --cell-size 10 --out triangular_mesh_1200.gts
```

and update `n_side` in a copy of `parameter1.txt` to 120 (1200/10). Compare
station time series between the two runs — the frictional patch itself
(Omega_f, the inner 800m x 800m) should be governed by rate-and-state
friction in the larger mesh too (only what's genuinely *outside* Omega_f
is locked at V=0), so this also exercises the `f_coef` friction-law
zoning; a version of `3dtri_BP8.f90` doing this would need an `n_side_omega_f`
parameter distinct from the mesh's own `n_side` to zone which elements are
frictional vs. locked — not yet implemented here, since the base case
(mesh = Omega_f exactly) has no locked elements to zone in the first
place.
