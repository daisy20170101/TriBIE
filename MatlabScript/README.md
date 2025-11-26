# MATLAB Triangular Green's Function Matrix Computation

This directory contains MATLAB scripts for computing the elastic Green's function (stiffness matrix) for triangular fault meshes using the Nikkhoo & Walter (2015) triangular dislocation method.

## Overview

The stiffness matrix represents the stress response at observation points (element centroids) due to unit slip on source triangular elements. This is used in boundary element method (BEM) simulations for earthquake cycle modeling.

## Files

- **`compute_trigreen_matrix.m`** - Main function for computing the stiffness matrix
- **`example_run.m`** - Example usage script with multiple configurations
- **`README.md`** - This file

## Requirements

- MATLAB (R2016b or later recommended)
- `TDstressHS.m` and `TDstressFS.m` from `../NikkhooWalter2015/` directory
- Triangular mesh file in GTS format

## Input Format

### GTS Mesh File Format

The input mesh file should be in GTS (GNU Triangulated Surface) format:

```
n_vertex n_edge n_cell
x1 y1 z1
x2 y2 z2
...
v1 v2 v3
v1 v2 v3
...
```

Where:
- Line 1: Number of vertices, edges, and cells
- Next `n_vertex` lines: Vertex coordinates (x, y, z in km)
- Next `n_cell` lines: Cell connectivity (vertex indices, 1-based)

### Example:

```
4 6 2
0.0 0.0 -5.0
1.0 0.0 -5.0
0.0 1.0 -5.0
1.0 1.0 -5.0
1 2 3
2 4 3
```

## Output Format

The script generates the following files:

### `trigreen_<ncore>.bin` Files

- One binary file per processor (numbered 0 to ncore-1)
- Contains `Nt × Ncell` double precision values
  - `Nt` = number of elements assigned to this processor
  - `Ncell` = total number of elements in the mesh
- Format: Row-major, each row contains influence coefficients for one source element on all observation elements
- Units: Bar (0.1 MPa)

### `position.bin` File

- Contains centroid positions of all elements
- Format: `Ncell × 3` double precision values (x, y, z coordinates)
- Written in column-major format (Fortran-compatible)

## Usage

### Basic Usage

```matlab
% Add path to Nikkhoo-Walter functions
addpath('../NikkhooWalter2015');

% Compute stiffness matrix with default parameters
compute_trigreen_matrix('triangular_mesh.gts', 4);
```

This will:
1. Read the mesh from `triangular_mesh.gts`
2. Distribute work across 4 processors
3. Generate files: `trigreen_0.bin`, `trigreen_1.bin`, `trigreen_2.bin`, `trigreen_3.bin`, and `position.bin`

### Advanced Usage with Custom Parameters

```matlab
compute_trigreen_matrix('mesh.gts', 8, ...
    'mu', 32000, ...          % Shear modulus (MPa)
    'nu', 0.28, ...           % Poisson's ratio
    'slip_ss', -1.0, ...      % Strike-slip component
    'slip_ds', 0.0, ...       % Dip-slip component
    'slip_ts', 0.0, ...       % Tensile-slip component
    'output_dir', './output'); % Output directory
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `mesh_file` | (required) | Path to GTS mesh file |
| `ncore` | (required) | Number of processors (load distribution) |
| `mu` | 30000 | Shear modulus (MPa) |
| `nu` | 0.25 | Poisson's ratio |
| `slip_ss` | -1.0 | Strike-slip component |
| `slip_ds` | 0.0 | Dip-slip component |
| `slip_ts` | 0.0 | Tensile-slip component |
| `output_dir` | `.` | Output directory path |

### Slip Convention

- **Strike-slip (ss)**: Positive for left-lateral, negative for right-lateral
- **Dip-slip (ds)**: Positive for reverse (thrust), negative for normal
- **Tensile-slip (ts)**: Positive for opening, negative for closing

## Load Balancing

The script automatically distributes elements across processors:

- **Even distribution**: If `n_cell` is divisible by `ncore`, each processor gets `n_cell/ncore` elements
- **Uneven distribution**: Extra elements are distributed to the first processors
  - Example: 100 cells, 8 processors → Processors 0-3 get 13 cells, processors 4-7 get 12 cells

## Performance Considerations

### Computational Complexity

- Total computations: `Ncell × Ncell` stress calculations
- Complexity: O(N²) where N = number of elements
- Each stress calculation involves solving the triangular dislocation problem

### Memory Requirements

For a mesh with `Ncell` elements and `ncore` processors:
- Memory per processor: ~ `(Ncell/ncore) × Ncell × 8 bytes`
- Example: 10,000 elements, 8 processors → ~100 MB per processor

### Computational Time

Approximate timing (depends on hardware):
- 100 elements: ~1 minute
- 1,000 elements: ~1-2 hours
- 10,000 elements: ~4-5 days (recommend parallel processing)

**Tip**: Use `ncore > 1` to distribute computation, even on a single machine (serial processing of chunks)

## Example Workflow

```matlab
% 1. Add required paths
addpath('../NikkhooWalter2015');

% 2. Define mesh and parameters
mesh_file = 'my_fault.gts';
ncore = 4;

% 3. Set material properties (example for crustal rocks)
mu = 32000;      % 32 GPa shear modulus
nu = 0.28;       % Poisson's ratio 0.28

% 4. Compute stiffness matrix
compute_trigreen_matrix(mesh_file, ncore, ...
    'mu', mu, ...
    'nu', nu, ...
    'output_dir', './trigreen_output');

% 5. Verify output
fprintf('Output files created:\n');
for i = 0:ncore-1
    fprintf('  trigreen_%d.bin\n', i);
end
fprintf('  position.bin\n');
```

## Comparison with Fortran Code

This MATLAB implementation produces the same output format as the Fortran code in `../TriGreen/calc_trigreen.f90`:

| Feature | Fortran (MPI+OpenMP) | MATLAB |
|---------|---------------------|--------|
| Parallelization | True parallel (MPI) | Sequential (load distribution) |
| Output format | Binary (stream) | Binary (compatible) |
| File naming | `trigreen_<myid>.bin` | `trigreen_<ncore>.bin` |
| Performance | Fast (compiled) | Slower (interpreted) |
| Ease of use | Requires compilation | Direct execution |

**Use cases**:
- **MATLAB**: Prototyping, small meshes (<1000 elements), development, verification
- **Fortran**: Production runs, large meshes (>1000 elements), HPC clusters

## Troubleshooting

### Common Issues

1. **"TDstressHS not found"**
   - Solution: Add path to NikkhooWalter2015 directory: `addpath('../NikkhooWalter2015')`

2. **"Cannot open mesh file"**
   - Solution: Check file path is correct and file exists

3. **"Out of memory"**
   - Solution: Increase `ncore` to reduce memory per processor, or use Fortran version

4. **"Results contain NaN"**
   - Possible causes: Observation points on triangle edges/vertices, degenerate triangles
   - Solution: Check mesh quality, ensure no duplicate vertices or zero-area triangles

### Verification

To verify output correctness:

```matlab
% Load output from one processor
fid = fopen('trigreen_0.bin', 'r');
data = fread(fid, [n_cell, Inf], 'double');
fclose(fid);

% Check for NaN or Inf
fprintf('NaN count: %d\n', sum(isnan(data(:))));
fprintf('Inf count: %d\n', sum(isinf(data(:))));
fprintf('Min value: %.6e\n', min(data(:)));
fprintf('Max value: %.6e\n', max(data(:)));
```

## References

1. Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical, artefact-free solution. *Geophysical Journal International*, 201(2), 1119-1141.

2. Stuart W.D., 1974. Angular dislocation in a half space. *Bulletin of the Seismological Society of America*, 64(6), 1851-1853.

## Support

For questions or issues:
- Check the main TriBIE documentation
- Review `example_run.m` for usage examples
- Consult Nikkhoo & Walter (2015) paper for method details

## License

This code is part of the TriBIE (Triangular Boundary Integral Element) package.
See main repository for license information.

---

**Author**: TriBIE Development Team
**Date**: 2025-11-18
**Version**: 1.0
