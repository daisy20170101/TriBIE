# Nikkhoo & Walter (2015) Triangular Dislocation Method

Fortran 90 implementation of the Nikkhoo & Walter (2015) method for calculating stresses and strains associated with triangular dislocations in an elastic half-space, with MPI/OpenMP parallel support for stiffness matrix calculation.

## Directory Structure

```
NikkhooWalter2015/
├── sub_nikkhoo.f90          # Core module with triangular dislocation calculations
├── calc_nikkhoo.f90         # MPI/OpenMP parallel stiffness calculation program
├── m_nikkhoo_green.f90      # Parameter module for parallel calculation
├── Makefile                 # Makefile for parallel build
├── Makefile.nikkhoo         # Makefile for serial test builds
├── runcompile_nikkhoo.sh    # Compilation script for parallel program
├── README_NIKKHOO.md        # This documentation
└── debug/                   # Debug and testing scripts
    ├── debug_nikkhoo.f90    # Fortran debug program
    ├── debug_matlab.m       # MATLAB reference for comparison
    ├── test_nikkhoo.f90     # Test program
    ├── Makefile.debug       # Debug build makefile
    └── compile_test.sh      # Compilation test script
```

## Reference

Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical, artefact-free solution. Geophysical Journal International.

## Features

- **Half-space solution**: Calculates stresses and strains in an elastic half-space
- **Triangular dislocations**: Supports arbitrary triangular fault elements
- **Multiple slip components**: Strike-slip, dip-slip, and tensile-slip
- **Artefact-free**: Uses the improved method to avoid numerical artifacts
- **MPI/OpenMP hybrid parallelization**: Efficient stiffness matrix calculation
- **Validated**: Results match MATLAB implementation (minor floating-point differences expected)

## Parallel Stiffness Calculation (calc_nikkhoo)

The `calc_nikkhoo` program calculates the Green's function stiffness matrix for triangular mesh elements using MPI for distributed computing and OpenMP for shared-memory parallelization within each MPI process.

### Quick Start

```bash
# Build
make

# Run with single process
./calc_nikkhoo

# Run with MPI (4 processes)
mpirun -np 4 ./calc_nikkhoo

# Run with hybrid MPI+OpenMP (4 MPI processes x 4 OpenMP threads each)
OMP_NUM_THREADS=4 mpirun -np 4 ./calc_nikkhoo
```

### Input Files

- `triangular_mesh.gts` - GTS format mesh file with triangular elements

### Output Files

- `nikkhoo_<rank>.bin` - Binary stiffness matrix file for each MPI rank
- `position.bin` - Centroid positions of all triangular elements

### Compilation Options

Using Make:
```bash
make                  # Build executable
make clean            # Remove build files
make cleanall         # Remove build and output files
make run              # Build and run
make run-mpi          # Build and run with 4 MPI processes
make run-mpi-8        # Build and run with 8 MPI processes
make help             # Show all options
```

Using shell script:
```bash
./runcompile_nikkhoo.sh
```

### Load Balancing

The program uses dynamic load balancing to distribute cells evenly across MPI processes:
- Base cells per process: `n_cell / size`
- Extra cells distributed to first `mod(n_cell, size)` processes
- Each process computes stiffness for its local cells against all source cells

## Core Module Functions (sub_nikkhoo.f90)

| Function | Description |
|----------|-------------|
| `tdstress_hs` | Main half-space stress/strain calculation |
| `tdstress_fs` | Full-space stress/strain calculation |
| `tdstress_harfunc` | Harmonic function contribution |
| `trimode_finder` | Point-in-triangle configuration detection |
| `tdsetup_s` | Angular dislocation setup and strain calculation |
| `angdis_strain` | Angular dislocation strain calculation |
| `angsetup_fsc_s` | Angular dislocation setup for free surface correction |
| `tens_trans` | Tensor transformation |
| `coord_trans` | Coordinate transformation between reference frames |

## Usage (Direct Module)

### Basic Usage

```fortran
use nikkhoo_walter
implicit none

real(DP) :: x, y, z
real(DP), dimension(3) :: p1, p2, p3
real(DP) :: ss, ds, ts, mu, lambda
real(DP), dimension(6) :: stress, strain

! Set up calculation point (must have z < 0)
x = 7.0_DP
y = -1.0_DP
z = -5.0_DP

! Set up triangular dislocation vertices
p1 = [-1.0_DP, -1.0_DP, -5.0_DP]
p2 = [1.0_DP, -1.0_DP, -5.0_DP]
p3 = [-1.0_DP, 1.0_DP, -4.0_DP]

! Set slip components
ss = 1.0_DP   ! Strike-slip
ds = -1.0_DP  ! Dip-slip
ts = 2.0_DP   ! Tensile-slip

! Set elastic parameters
mu = 3.0e10_DP      ! Shear modulus
lambda = 3.0e10_DP  ! Lame's first parameter

! Calculate stresses and strains
call tdstress_hs(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)
```

### Output Format

**Stress tensor** (6 components):
- `stress(1:6)` = Sxx, Syy, Szz, Sxy, Sxz, Syz

**Strain tensor** (6 components):
- `strain(1:6)` = Exx, Eyy, Ezz, Exy, Exz, Eyz

## Configuration Detection

The `trimode_finder` function determines the point's position relative to the triangle:

| trimode | Configuration | Description |
|---------|--------------|-------------|
| 1 | I | Point inside triangle projection |
| -1 | II | Point outside triangle projection |
| 0 | III | Point on triangle edge (returns NaN) |

**Note**: Points exactly on triangle edges (trimode=0) return NaN for strain/stress values. This matches the MATLAB implementation and is expected behavior for singular points.

## Important Notes

1. **Half-space constraint**: All Z coordinates must be negative (below the free surface)
2. **Coordinate system**: Uses East-North-Up (ENU) coordinate system
3. **Units**: Input coordinates and elastic parameters should be in consistent units
4. **Edge points**: Points on triangle edges return NaN (singular)
5. **Numerical precision**: Uses double precision; minor floating-point differences from MATLAB are expected

## Material Parameters (m_nikkhoo_green.f90)

Default parameters can be modified in `m_nikkhoo_green.f90`:

```fortran
! Subduction geometry
real(DP), parameter :: subd_az = 0.d0   ! Subduction azimuth (degree)
real(DP), parameter :: rot_deg = 0.d0  ! Rotation angle

! Material parameters
real(DP), parameter :: parm_nu = ...    ! Poisson's ratio
real(DP), parameter :: parm_miu = 30000.d0  ! Rigidity (MPa)
real(DP), parameter :: parm_l = ...     ! Lambda (Lame parameter, MPa)
```

## Testing and Validation

### Run Debug Program

```bash
cd debug
gfortran -o debug_nikkhoo debug_nikkhoo.f90 ../sub_nikkhoo.f90
./debug_nikkhoo
```

### Compare with MATLAB

```matlab
% Run debug_matlab.m in MATLAB and compare with Fortran output
run('debug/debug_matlab.m')
```

## Recent Updates (December 2024)

### Bug Fixes

1. **trimode_finder Index Mapping**
   - Issue: Incorrect array index mapping for barycentric coordinates
   - Fix: Changed `p(2)` to `p(3)` for z-coordinates (MATLAB uses 2-element vectors, Fortran uses 3-element)

2. **Transformation Matrix Construction**
   - Issue: Inconsistent matrix construction between functions
   - Fix: Standardized to column-based construction (`A(:,i) = vec`) for use with `coord_trans`

3. **Edge Point Handling (casez_log)**
   - Issue: Attempted to compute averages for singular edge points
   - Fix: Now returns NaN for edge points, matching MATLAB `TDstressFS.m` behavior

### New Features

- **MPI/OpenMP Parallel Stiffness Calculation** (`calc_nikkhoo.f90`)
  - Hybrid parallelization for large-scale problems
  - Dynamic load balancing across MPI processes
  - OpenMP parallelization within each process
  - Compatible with TriGreen output format

## Comparison with TriGreen

| Feature | TriGreen (Stuart) | NikkhooWalter2015 |
|---------|-------------------|-------------------|
| Method | Stuart (1988) | Nikkhoo & Walter (2015) |
| Domain | Full-space approximation | True half-space |
| Artifacts | Potential near-field artifacts | Artefact-free |
| Main program | `calc_trigreen.f90` | `calc_nikkhoo.f90` |
| Output | `trigreen_<rank>.bin` | `nikkhoo_<rank>.bin` |

## Integration with TriBIE

The stiffness matrices generated by `calc_nikkhoo` are compatible with the main TriBIE simulation framework. Output files follow the same format as `calc_trigreen` for seamless integration.
