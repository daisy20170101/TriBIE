# Python Implementation of Triangular Dislocation Stress Calculations

Python translation of the MATLAB code for triangular dislocation stress/strain calculations in elastic media.

## Reference

**Nikkhoo M. and Walter T.R., 2015.** Triangular dislocation: An analytical, artefact-free solution. *Geophysical Journal International*.

Original MATLAB code: Mehdi Nikkhoo (mehdi.nikkhoo@gmail.com)
Python translation: From Fortran debugging session (2025)

## Features

- **Full-Space Solution** (`tdstress_fs`): Calculate stress and strain for triangular dislocations in infinite elastic medium
- **Half-Space Solution** (`tdstress_hs`): Calculate stress and strain accounting for free surface (Note: harmonic function contribution simplified)
- **Numpy-based**: Efficient array operations matching MATLAB behavior
- **Well-documented**: Docstrings for all functions with parameter descriptions

## Installation

```python
# Add the python_tdstress directory to your Python path
import sys
sys.path.insert(0, '/path/to/TriBIE/python_tdstress')

# Import the main functions
from tdstress_fs import tdstress_fs
from tdstress_hs import tdstress_hs
```

Or install as a package:
```bash
cd /path/to/TriBIE/python_tdstress
pip install -e .
```

## Usage

### Basic Example - Full Space

```python
import numpy as np
from tdstress_fs import tdstress_fs

# Define triangle vertices (East, North, Up coordinates)
P1 = np.array([-1.0, -1.0, -5.0])
P2 = np.array([1.0, -1.0, -5.0])
P3 = np.array([-1.0, 1.0, -4.0])

# Slip components
Ss = 1.0   # Strike-slip
Ds = -1.0  # Dip-slip
Ts = 2.0   # Tensile-slip

# Elastic parameters (Lame constants)
mu = 3.0e10     # Shear modulus
lam = 3.0e10    # Lame's first parameter

# Calculation point
X, Y, Z = 0.0, 0.0, -5.0

# Calculate stress and strain
stress, strain = tdstress_fs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, mu, lam)

# Access components
print(f"Exx = {strain[0, 0]:.6e}")
print(f"Stress tensor: {stress[0]}")
```

### Multiple Points

```python
# Calculate for multiple points
X = np.array([0.0, 1.0, 2.0])
Y = np.array([0.0, 0.0, 0.0])
Z = np.array([-5.0, -5.0, -5.0])

stress, strain = tdstress_fs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, mu, lam)

# Results are (n_points, 6) arrays
print(f"Shape: {strain.shape}")  # (3, 6)
```

### Half-Space Solution

```python
from tdstress_hs import tdstress_hs

# All Z coordinates must be negative (below surface)
stress_hs, strain_hs = tdstress_hs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, mu, lam)
```

## Coordinate Systems

### EFCS (Earth-Fixed Coordinate System)
- **X**: East
- **Y**: North
- **Z**: Up (must be negative for half-space)

### TDCS (Triangular Dislocation Coordinate System)
- **Origin**: At vertex P2
- **x-axis**: Normal to triangle (perpendicular)
- **y-axis, z-axis**: In the plane of the triangle

Transformations between coordinate systems are handled automatically.

## Output Format

Both `tdstress_fs` and `tdstress_hs` return:

**Stress** : ndarray, shape (n_points, 6)
- Components: [Sxx, Syy, Szz, Sxy, Sxz, Syz]
- Units: Same as input Lame constants

**Strain** : ndarray, shape (n_points, 6)
- Components: [Exx, Eyy, Ezz, Exy, Exz, Eyz]
- Units: Dimensionless

## Modules

### `td_utils.py`
Helper functions for coordinate and tensor transformations:
- `coord_trans`: Transform coordinates between systems
- `tens_trans`: Transform tensor components
- `trimodefinder`: Determine configuration mode (inside/outside/on-edge)

### `ang_dislocation.py`
Angular dislocation calculations:
- `ang_dis_strain`: Calculate strain for angular dislocation
- `td_setup_s`: Transform and calculate in local coordinate system

### `tdstress_fs.py`
Full-space triangular dislocation implementation

### `tdstress_hs.py`
Half-space triangular dislocation implementation
- Includes main dislocation + image dislocation
- **Note**: Harmonic function contribution simplified (returns zeros)

## Testing

Run the test script to verify against reference values:

```bash
cd python_tdstress
python test_tdstress.py
```

Expected output:
```
Testing Triangular Dislocation Implementation
======================================================================
...
Point 1 (center): (-0.333, -0.333, -4.667)
  Expected Exx: 4.810470052551810e-02
  Got Exx:      4.810470052551810e-02
  Error:        0.000000000000000e+00
  Rel. Error:   0.000000%
  *** PASS ***
...
```

## Known Limitations

1. **Harmonic Function**: The half-space solution (`tdstress_hs`) currently has a simplified implementation of the harmonic function contribution. It returns zeros instead of the full free-surface correction. This means results will match `TDstressFS` + image dislocation but not the complete `TDstressHS` solution.

2. **Singular Points**: Points on triangle edges in the triangle plane return NaN (mathematically correct but may need regularization for applications).

3. **Performance**: No special optimizations for large arrays yet (straightforward NumPy translation).

## Future Improvements

- [ ] Complete implementation of harmonic function contribution (`AngSetupFSC_S`)
- [ ] Vectorized operations for better performance
- [ ] Optional regularization for near-singular points
- [ ] Additional validation tests against MATLAB/Fortran
- [ ] Proper Python package setup with `setup.py`

## File Structure

```
python_tdstress/
├── __init__.py              # Package initialization
├── td_utils.py              # Coordinate/tensor transformations, trimode
├── ang_dislocation.py       # Angular dislocation calculations
├── tdstress_fs.py          # Full-space solution
├── tdstress_hs.py          # Half-space solution
├── test_tdstress.py        # Test script
└── README.md               # This file
```

## License

Copyright (c) 2014 Mehdi Nikkhoo (original MATLAB code)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Contact

For questions about the original MATLAB implementation:
- Mehdi Nikkhoo: mehdi.nikkhoo@gmail.com

For questions about this Python translation:
- See TriBIE repository issues
