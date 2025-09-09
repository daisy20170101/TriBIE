# Nikkhoo & Walter (2015) Triangular Dislocation Method

This directory contains a Fortran 90 implementation of the Nikkhoo & Walter (2015) method for calculating stresses and strains associated with triangular dislocations in an elastic half-space.

## Files

- `sub_nikkhoo.f90` - Main module containing the triangular dislocation calculations
- `test_nikkhoo.f90` - Test program demonstrating usage
- `Makefile.nikkhoo` - Makefile for compiling the test program
- `README_NIKKHOO.md` - This documentation file

## Reference

Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical, artefact-free solution. Geophysical Journal International.

## Features

- **Half-space solution**: Calculates stresses and strains in an elastic half-space
- **Triangular dislocations**: Supports arbitrary triangular fault elements
- **Multiple slip components**: Strike-slip, dip-slip, and tensile-slip
- **Artefact-free**: Uses the improved method to avoid numerical artifacts
- **Fortran 90**: Modern Fortran implementation with modules

## Usage

### Basic Usage

```fortran
use nikkhoo_walter
implicit none

integer, parameter :: n_points = 100
real(DP), dimension(n_points) :: x, y, z
real(DP), dimension(3) :: p1, p2, p3
real(DP) :: ss, ds, ts, mu, lambda
real(DP), dimension(n_points, 6) :: stress, strain

! Set up calculation points
x = [your x coordinates]
y = [your y coordinates]  
z = [your z coordinates]

! Set up triangular dislocation vertices
p1 = [x1, y1, z1]
p2 = [x2, y2, z2]
p3 = [x3, y3, z3]

! Set slip components
ss = 1.0_DP  ! Strike-slip
ds = 0.5_DP  ! Dip-slip
ts = 0.0_DP  ! Tensile-slip

! Set elastic parameters
mu = 3.0e10_DP      ! Shear modulus
lambda = 3.0e10_DP  ! Lame's first parameter

! Calculate stresses and strains
call tdstress_hs(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain, n_points)
```

### Output Format

The output arrays contain:
- `stress(:, 1)` - Sxx component
- `stress(:, 2)` - Syy component  
- `stress(:, 3)` - Szz component
- `stress(:, 4)` - Sxy component
- `stress(:, 5)` - Sxz component
- `stress(:, 6)` - Syz component

- `strain(:, 1)` - Exx component
- `strain(:, 2)` - Eyy component
- `strain(:, 3)` - Ezz component
- `strain(:, 4)` - Exy component
- `strain(:, 5)` - Exz component
- `strain(:, 6)` - Eyz component

## Compilation

### Using the provided Makefile

```bash
# Compile the test program
make -f Makefile.nikkhoo

# Run the test
make -f Makefile.nikkhoo run

# Clean up
make -f Makefile.nikkhoo clean
```

### Manual compilation

```bash
# Compile the module
gfortran -c sub_nikkhoo.f90

# Compile your program
gfortran -c your_program.f90

# Link everything
gfortran -o your_program your_program.o sub_nikkhoo.o
```

## Integration with calc_trigreen

To integrate this with your existing `calc_trigreen` code:

1. **Replace the stiffness calculation**: Use `tdstress_hs` instead of the current `sub_comdun.f90` function
2. **Modify the interface**: Adapt the function calls to match your existing data structures
3. **Update the Makefile**: Add `sub_nikkhoo.f90` to your compilation

### Example integration

```fortran
! In your main calc_trigreen code
use nikkhoo_walter

! Replace your existing stiffness calculation with:
call tdstress_hs(x_points, y_points, z_points, &
                 p1, p2, p3, &
                 slip_strike, slip_dip, slip_tensile, &
                 shear_modulus, lame_parameter, &
                 stress_tensor, strain_tensor, n_calc_points)
```

## Important Notes

1. **Half-space constraint**: All Z coordinates must be negative (below the free surface)
2. **Coordinate system**: Uses East-North-Up (ENU) coordinate system
3. **Units**: Input coordinates and elastic parameters should be in consistent units
4. **Performance**: The current implementation includes placeholder functions for some complex calculations

## Limitations

- Some functions (`angdis_strain`, `angsetup_fsc_s`) are currently placeholder implementations
- The full MATLAB functionality has not been completely ported yet
- Performance optimizations may be needed for large-scale calculations

## Future Improvements

- Complete implementation of all strain calculation functions
- Add OpenMP parallelization for better performance
- Optimize memory usage for large calculations
- Add more comprehensive error checking and validation

## Testing

Run the test program to verify the installation:

```bash
make -f Makefile.nikkhoo run
```

This will calculate stresses and strains for a simple triangular dislocation and display the results.

## Support

For questions or issues with this implementation, please refer to the original MATLAB code or the Nikkhoo & Walter (2015) paper for reference.
