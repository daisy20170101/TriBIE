# Nikkhoo & Walter (2015) Triangular Dislocation Method

Fortran 90 implementation of the Nikkhoo & Walter (2015) method for calculating stresses and strains associated with triangular dislocations in an elastic half-space.

## Directory Structure

```
NikkhooWalter2015/
├── sub_nikkhoo.f90      # Main module with all triangular dislocation calculations
├── Makefile.nikkhoo     # Makefile for building programs
├── README_NIKKHOO.md    # This documentation
└── debug/               # Debug and testing scripts
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
- **Validated**: Results match MATLAB implementation (minor floating-point differences expected)

## Core Functions

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

## Usage

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

## Compilation

```bash
# Compile the module
gfortran -c sub_nikkhoo.f90

# Link with your program
gfortran -o your_program your_program.f90 sub_nikkhoo.o
```

### Using the Makefile

```bash
make -f Makefile.nikkhoo
```

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

## Recent Bug Fixes (December 2024)

### 1. trimode_finder Index Mapping
- **Issue**: Incorrect array index mapping for barycentric coordinates
- **Fix**: Changed `p(2)` to `p(3)` for z-coordinates (MATLAB uses 2-element vectors, Fortran uses 3-element)

### 2. Transformation Matrix Construction
- **Issue**: Inconsistent matrix construction between functions
- **Fix**: Standardized to column-based construction (`A(:,i) = vec`) for use with `coord_trans`

### 3. Edge Point Handling (casez_log)
- **Issue**: Attempted to compute averages for singular edge points
- **Fix**: Now returns NaN for edge points, matching MATLAB `TDstressFS.m` behavior

## Integration with TriBIE

To integrate with the main TriBIE stiffness calculation:

```fortran
use nikkhoo_walter

! Replace existing stiffness calculation:
do i = 1, n_calc_points
    call tdstress_hs(x_points(i), y_points(i), z_points(i), &
                     p1, p2, p3, &
                     slip_strike, slip_dip, slip_tensile, &
                     shear_modulus, lame_parameter, &
                     stress_tensor(i, :), strain_tensor(i, :))
end do
```
