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

## Mesh: any size >= Omega_f; the true frictional domain is fixed

`3dtri_BP8.f90` zones each element by its centroid: inside Omega_f
(|x2|<400m and |x3|<400m, fixed by BP8 Table 1 -- `lf_fixed` in
`module_bp8.f90`, a compile-time constant, *not* read from
parameter1.txt) it gets rate-and-state friction as normal; outside it,
Eq. (13)'s V=0 boundary condition is enforced exactly (its state is never
integrated). The mesh itself just needs to cover at least Omega_f --
`example4/stiffness/triangular_mesh.gts` here is the default 800m x 800m
case (mesh = exactly Omega_f, no extra elements), but a bigger mesh
(e.g. `make_bp8_mesh.py --domain 1200`) works too, with the extra
elements automatically locked. `read_parameters` errors out if the mesh
is smaller than Omega_f.

Locked elements still matter for output: BP8 Eq. (8) says their shear
stress is `tau0 + Dtau(t)` (V=0 removes the radiation-damping term, not
the elastic-coupling one) -- it is *not* simply zero, and reporting zero
there corrupts bilinear-interpolated output (profile lines, stations) near
the Omega_f boundary on any mesh bigger than Omega_f. `derivs()` tracks
`Dtau2`/`Dtau3` for locked elements by repurposing their (otherwise
unused, since V=0) velocity state slots; `write_all_output` reports
`tau0 + Dtau` for them instead of zero. This was a real bug caught by
diffing a padded (1000m) mesh against an exact (800m) one at the same
resolution: shear stress near the boundary was off by up to 50% before
the fix, and matches to ~1e-6 relative (floating-point-level) after it.

Profile-line output (Section 4.3) is always the fixed 81-node, exactly-
10m-spaced grid from -400 to 400 regardless of the mesh's own resolution
or domain size (`n_nodes_fixed`/`node_dz_fixed` in `module_bp8.f90`) --
interpolated via `bilinear_cell` from whatever the mesh actually provides.

The benchmark's own tip (Section 6) is to verify results are independent
of computational domain size; see "Domain-independence check" below for
how to actually exercise the padded-mesh path now that it's implemented.

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
| 3 | `n_side dz_cell` — mesh cells across the full domain, cell size (m). Must satisfy `n_side*dz_cell >= 800` (mesh covers at least Omega_f, lf=400m fixed) or the run errors out at startup. |
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

To verify that meshing beyond Omega_f doesn't change the result, rerun
steps 1-3 with a larger domain, e.g.:

```bash
mkdir -p stiffness_1200 && cd stiffness_1200
python3 ../make_bp8_mesh.py --domain 1200 --cell-size 10 --out triangular_mesh.gts
mpirun -np 4 ../../NikkhooWalter2015/calc_nikkhoo_fs
```

then a copy of `parameter1.txt` with `n_side=120` (1200/10) and
`stiffname=./stiffness_1200/`. No other changes are needed — the extra
elements outside |x2|,|x3|<400 are automatically locked (Eq. 13), and
profile/station output still reports the fixed -400..400, 10m-spaced grid
either way.

This was actually exercised during development: a 1000m/10m mesh (n_side=100)
compared against the exact 800m/10m case matched to ~1e-6 relative
(floating-point level) in shear stress, slip, and pore pressure profiles,
and `global.dat` (Vmax, moment rate) matched exactly. If you rerun this
check yourself and see a bigger discrepancy than that, look at
`is_active`/`tau0_2`/`tau0_3` in `3dtri_BP8.f90` first — that's where the
padded-mesh-specific logic lives.
