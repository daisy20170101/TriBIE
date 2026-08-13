# DynamicRuptureInput - TriBIE Input Files Repository

## Overview
This repository contains all the input files required for running dynamic rupture simulations using the TriBIE (Triangular Boundary Integral Element) code. TriBIE is a Fortran90 parallel (MPI) code designed to simulate slow slip events and earthquake cycles in arbitrary curved faults buried in half-space media.

## Repository Structure
```
DynamicRuptureInput/
├── README.md                    # This file
├── parameter1.txt              # Main simulation configuration
├── var-BP5_h500_140_60.dat    # On-fault physical parameters
├── profdp-BP5_h500_140_60.dat # Dip direction observation points
├── profstrk-BP5_h500_140_60.dat # Strike direction observation points
├── obvs.dat                    # Observation points configuration
├── input_cal_stiffness.txt     # Stiffness calculation parameters
├── fault_h500_140_60.gts      # Triangular fault mesh (GTS format)
└── area-BP5_h500_140_60.dat   # Element areas for moment calculation
```

## File Descriptions

### 1. Main Configuration File

#### `parameter1.txt`
The primary configuration file that controls all aspects of the simulation.

**Format:**
```
<jobname>                    ! Line 1: Simulation job identifier
<foldername>                 ! Line 2: Output directory path
<stiffname>                  ! Line 3: Stiffness matrix file prefix
<restartname>                ! Line 4: Restart file name (if applicable)
<Nab> <Nt_all> <Nt> <Lratio> <nprocs> <n_obv> <np1> <np2>  ! Line 5: Array dimensions
<Idin> <Idout> <Iprofile> <Iperb> <Isnapshot>               ! Line 6: Control flags
<Vpl>                        ! Line 7: Plate velocity (m/s)
<tmax>                       ! Line 8: Maximum simulation time (years)
<tslip_ave> <tslipend> <tslip_aveint>                       ! Line 9: Slip averaging parameters
<tint_out> <tmin_out> <tint_cos> <tint_sse>                 ! Line 10: Output intervals
<vcos> <vsse1> <vsse2>                                      ! Line 11: Velocity thresholds
<nmv> <nas> <ncos> <nnul> <nsse> <n_nul_int>               ! Line 12: Output counters
<s1(1)> <s1(2)> ... <s1(10)>                               ! Line 13: Additional parameters
```

**Parameter Details:**
- **`jobname`**: Unique identifier for the simulation run
- **`foldername`**: Directory where output files will be saved
- **`stiffname`**: Prefix for stiffness matrix files (e.g., "trigreen_")
- **`restartname`**: Name of restart file for continuing previous simulation
- **`Nab`**: Number of temperature-dependent a-b profiles
- **`Nt_all`**: Total number of fault elements
- **`Nt`**: Local number of elements per process (auto-calculated)
- **`nprocs`**: Number of MPI processes
- **`Vpl`**: Plate velocity (typically 1e-10 to 1e-9 m/s)
- **`tmax`**: Maximum simulation time in years

### 2. Physical Parameters

#### `var-BP5_h500_140_60.dat`
Contains the physical properties for each fault element.

**Format:**
```
Column 1: Effective normal stress (Pa)
Column 2: Characteristic slip distance, dc (m)
Column 3: Rate-and-state friction parameter (a-b)
Column 4: Initial slip rate (m/s)
Column 5: Initial state variable (s)
```

**Notes:**
- One row per fault element
- Must match the number of elements in the mesh
- Values should be physically reasonable for your fault model

### 3. Observation Points

#### `profdp-BP5_h500_140_60.dat` (Dip Direction)
#### `profstrk-BP5_h500_140_60.dat` (Strike Direction)
Contains coordinates of observation points for monitoring fault behavior.

**Format:**
```
Column 1: X coordinate (m)
Column 2: Y coordinate (m) 
Column 3: Z coordinate (m)
```

**Usage:**
- Used for generating displacement and velocity profiles
- Coordinates should be in the same reference system as the fault mesh
- Can be placed on or off the fault surface

#### `obvs.dat`
Configuration file for observation points.

**Format:**
```
Line 1: Number of observation points
Subsequent lines: Observation point configurations
```

### 4. Mesh Files

#### `fault_h500_140_60.gts`
Triangular fault mesh in GTS (GNU Triangulated Surface) format.

**Format:**
```
Line 1: <number_of_vertices> <number_of_edges> <number_of_cells>
Lines 2 to number_of_vertices+1: Vertex coordinates (x, y, z)
Remaining lines: Triangle topology (3 vertex indices per triangle)
```

**Notes:**
- GTS is a standard format for triangular surface meshes
- Can be generated using various mesh generation tools
- Element size should be appropriate for the physics being modeled

#### `area-BP5_h500_140_60.dat`
Element areas for moment calculation.

**Format:**
```
Column 1: Element area (m²)
```

**Usage:**
- Used for calculating seismic moment
- Must correspond to the mesh elements
- Areas are used in stress calculations

### 5. Stiffness Calculation

#### `input_cal_stiffness.txt`
Parameters for stiffness matrix calculation.

**Format:**
```
Line 1: Number of processes
Line 2: Mesh file path
Line 3: Output directory
```

## Usage Instructions

### 1. Prepare Input Files
1. **Mesh Generation**: Create or obtain a triangular fault mesh in GTS format
2. **Parameter Setup**: Configure `parameter1.txt` with your simulation parameters
3. **Physical Properties**: Set up `var-*.dat` with appropriate friction parameters
4. **Observation Points**: Define observation locations in profile files

### 2. Run Stiffness Calculation
```bash
cd TriGreen/
./runcompile.sh
mpirun -np <n_processes> ./calc_trigreen
```

### 3. Run Simulation
```bash
cd src/
./compile.sh
cd ../input
mpirun -np <n_processes> ../src/3dtri_BP5
```

## Example: SEAS BP5 Benchmark

The repository includes input files for the SEAS (Sequences of Earthquakes and Aseismic Slip) BP5 benchmark problem:

- **Fault Geometry**: 140 km × 60 km rectangular fault
- **Element Size**: 500 m triangular elements
- **Physics**: Rate-and-state friction with temperature-dependent parameters
- **Validation**: Verified against SCEC community standards

## File Naming Convention

Files follow the pattern: `{type}-{model}_{resolution}_{strike}_{dip}.{extension}`

- **`type`**: File purpose (var, prof, area, fault)
- **`model`**: Model identifier (BP5)
- **`resolution`**: Element size (h500 = 500m)
- **`strike`**: Strike length in km (140)
- **`dip`**: Dip length in km (60)

## Troubleshooting

### Common Issues
1. **File Not Found**: Ensure all input files are in the correct directory
2. **Parameter Mismatch**: Verify array dimensions match mesh element count
3. **MPI Errors**: Check that process count matches stiffness file count
4. **Mesh Issues**: Validate GTS file format and element connectivity

### Validation
- Check that `Nt_all` matches the number of elements in your mesh
- Verify that `nprocs` matches the number of stiffness files
- Ensure observation point coordinates are within reasonable bounds

## References

- **TriBIE Code**: [GitHub Repository](https://github.com/daisy20170101/TriBIE)
- **SEAS Project**: [SCEC Website](https://strike.scec.org/cvws/seas/)
- **GTS Format**: [GNU Triangulated Surface](http://gts.sourceforge.net/)

## Contact

For questions about these input files or the TriBIE code:
- **Author**: D. Li (d.li@gns.cri.nz)
- **Last Updated**: July 2024

## License

This repository is part of the TriBIE project. See the main repository for licensing information.
