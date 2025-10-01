# TriBIE

**TriBIE** (Triangular Boundary Integral Equation) is a high-performance Fortran90 parallel computing framework for simulating earthquake cycles, slow slip events, and aseismic transients on complex 3D fault geometries. The code employs hybrid MPI+OpenMP parallelization with advanced computational optimizations including SIMD vectorization and dynamic load balancing to efficiently model rate-and-state friction physics on triangular fault meshes embedded in elastic half-space media.

## Key Features

- **🌍 Complex Fault Geometries**: Supports arbitrary curved faults with triangular mesh discretization
- **⚡ High-Performance Computing**: Hybrid MPI+OpenMP parallelization with SIMD optimization
- **🔧 Advanced Physics**: Rate-and-state friction laws for realistic earthquake cycle modeling
- **📊 Modern I/O**: HDF5/XDMF output for direct visualization in Paraview
- **⚖️ Dynamic Load Balancing**: Automatic work distribution across irregular mesh geometries
- **🔄 Accumulative Simulations**: Long-term earthquake cycle studies with restart capability
- **✅ Scientifically Validated**: Verified in the SCEC Sequences of Earthquakes and Aseismic Slip (SEAS) Project

**Applications**: Earthquake cycle modeling, slow slip event analysis, tsunami hazard assessment, and fault system dynamics research.

**Verification**: TriBIE has been extensively validated through the [SCEC SEAS Project](https://strike.scec.org/cvws/seas/) benchmark comparisons.



## Simulation

![varBP5](https://github.com/daisy20170101/TriBIE/assets/33549997/c0b43d1b-777a-48e0-bda4-72c7a9b0e95e)
Figure: Mapview of on-fault distribution of fault key parameters in BP5 example.

## Results
![bp5_slip_h1000_140_60](https://github.com/daisy20170101/TriBIE/assets/33549997/b5e8804c-297d-4bbd-b8cf-c75a64f6bea9)
Figure: Cumulative slip along strike (left) and along downdip (right) in the first 800 modeling years. Coseismic slip in red while interseismic in blue. 

---

## User Guide

### Overview
TriBIE provides a complete workflow for earthquake cycle simulation in three main phases:
1. **Stiffness Matrix Computation** using `calc_trigreen.f90`
2. **Cycling Simulation** using `3dtri_BP5.f90` 
3. **Result Analysis** using output visualization tools

### Workflow Overview

```
Mesh Generation → Stiffness Computation → Cycling Simulation → Analysis
     ↓                    ↓                    ↓              ↓
  triangular_mesh.gts  calc_trigreen.f90  3dtri_BP5.f90   Results
```

### 1. Stiffness Matrix Computation

#### Prerequisites
- **Input Mesh**: `triangular_mesh.gts` file containing triangular fault elements
- **MPI Environment**: Multi-process execution environment
- **Compilation**: Use `TriGreen/runcompile.sh`

#### Execution
```bash
cd TriGreen/
./runcompile.sh
mpirun -np <n_processes> ./calc_trigreen
```

#### Output Files
- **`trigreen_<process_id>.bin`**: Stiffness matrix for each process
- **`position.bin`**: Element centroid positions
- **Console Output**: Distribution information and performance metrics

#### Key Features
- **Dynamic Load Balancing**: Automatically distributes work unevenly across processes
- **SIMD Optimization**: Vectorized computation for improved performance
- **Memory Management**: Optimized allocation/deallocation strategies

### 2. Cycling Simulation

#### Prerequisites
- **Stiffness Files**: `trigreen_<process_id>.bin` files from step 1
- **Parameter File**: `input/parameter1.txt` with simulation parameters
- **Compilation**: Use `src/compile.sh`

#### Execution
```bash
cd src/
./compile.sh
cd ../input
mpirun -np <n_processes> ../src/3dtri_BP5
```

#### Key Features
- **MPI_Scatterv**: Proper handling of uneven distributions
- **Rate-and-State Friction**: Physics-based fault behavior modeling
- **SIMD Optimization**: Vectorized physics calculations
- **Dynamic Load Balancing**: Compatible with calc_trigreen distribution

### 3. Parameter Configuration

#### `parameter1.txt` Format
```
<jobname>                    ! Simulation job identifier
<foldername>                 ! Output directory path
<stiffname>                  ! Stiffness matrix file prefix
<restartname>                ! Restart file name (if applicable)
<Nt_all> <nprocs> <n_obv> <num_of_receivers_along_strike> <num_of_receivers_along_downdip>  ! Array dimensions
<Idin> <Idout> <Iprofile> <Iperb> <Isnapshot>               ! Control flags
<Vpl>                        ! Plate velocity (m/s)
<tmax>                       ! Maximum simulation time (years)
<tslip_ave> <tslipend> <tslip_aveint>                       ! Slip averaging parameters
<tint_out> <tmin_out> <tint_cos> <tint_sse>                 ! Output intervals
<vcos> <vsse1> <vsse2>                                      ! Velocity thresholds
<nmv> <nas> <ncos> <nnul> <nsse> <n_nul_int>               ! Output counters
<s1(1)> <s1(2)> ... <s1(10)>                               ! Additional parameters
```

#### Key Parameter Descriptions

**Array Dimensions:**
- **`Nt_all`**: Total number of fault elements
- **`nprocs`**: Number of MPI processes
- **`n_obv`**: Number of observation points

**Physical Parameters:**
- **`Vpl`**: Plate velocity (typically 1e-10 to 1e-9 m/s)
- **`tmax`**: Maximum simulation time in years

**Output Control:**
- **`tint_cos`**: Coseismic output interval
- **`tint_sse`**: Slow slip event output interval

### 4. Example Workflow

#### Step 1: Prepare Mesh
```bash
# Ensure triangular_mesh.gts exists in TriGreen/ directory
ls TriGreen/triangular_mesh.gts
```

#### Step 2: Compute Stiffness
```bash
cd TriGreen/
./runcompile.sh
mpirun -np 8 ./calc_trigreen
# Generates: trigreen_0.bin, trigreen_1.bin, ..., trigreen_7.bin
```

#### Step 3: Configure Parameters
```bash
cd input/
# Edit parameter1.txt with your simulation parameters
# Ensure nprocs matches the number of stiffness files
```

#### Step 4: Run Simulation
```bash
cd src/
./compile.sh
cd ../input/
mpirun -np 8 ../src/3dtri_BP5
```

#### Step 5: Analyze Results
```bash
# Check output files for results
ls -la area* fltst* prof*
```

### 5. HDF5 and XDMF Output

TriBIE now supports modern HDF5/XDMF output for visualization:

#### Features
- **HDF5 Time Series**: Efficient storage of large temporal datasets
- **XDMF Visualization**: Direct compatibility with Paraview
- **Accumulative Writing**: Multiple simulation runs append to existing data
- **Mesh Integration**: Automatic mesh export for visualization

#### Output Files
- **`timeseries_data_<jobname>.h5`**: Main time series data
- **`timeseries_data_<jobname>.xdmf`**: Paraview visualization file
- **`sse_timeseries_data_<jobname>.h5`**: SSE-specific data
- **`sse_timeseries_data_<jobname>.xdmf`**: SSE visualization file

### 6. Troubleshooting

#### Common Issues

**MPI Communication Errors:**
- **Symptom**: `MPI_ERR_TRUNCATE` or communication failures
- **Solution**: Ensure `nprocs` in parameter1.txt matches the number of stiffness files

**Memory Issues:**
- **Symptom**: `free(): invalid size` or allocation failures  
- **Solution**: Proper memory management implemented with conditional deallocation

**HDF5 Errors:**
- **Symptom**: "name already exists" or dataset creation failures
- **Solution**: Fixed with existence checks for groups and datasets

#### Performance Optimization
- **SIMD**: Use `*_simd.f90` versions for better performance
- **OpenMP**: Hybrid MPI+OpenMP parallelization available
- **Process Count**: Match MPI processes to available cores

### 7. Advanced Features

#### Dynamic Load Balancing
- Automatically handles uneven element distributions
- Optimizes work distribution across processes
- Compatible with irregular fault geometries

#### SIMD Vectorization
- Automatic vectorization of computational loops
- Cache-aware memory access patterns
- Performance improvements on modern processors

#### Parallel I/O
- HDF5-based parallel file output
- XDMF integration for visualization
- Accumulative time series storage

### 8. Best Practices

#### Performance
- Use SIMD-optimized versions for production runs
- Match MPI processes to available hardware cores
- Monitor memory usage and adjust problem size accordingly
- Set appropriate OpenMP thread counts for hybrid parallelization

#### Reliability
- Always verify input file formats and parameters
- Use restart capability for long simulations
- Check output files for expected results
- Test with smaller problems before large-scale runs

For detailed parameter descriptions and advanced usage, see `src/USER_GUIDE.md`.

---

## Running Example1 - BP5 Benchmark

The `example1/` directory contains a complete working example based on the SCEC SEAS BP5 benchmark problem. This example demonstrates earthquake cycle simulation on a planar fault with rate-and-state friction.

### Example1 Overview

**Problem**: SCEC SEAS Benchmark Problem 5 (BP5) - Long-term earthquake cycles on a vertical strike-slip fault
- **Fault geometry**: 160 km × 60 km planar fault 
- **Depth**: Surface to 60 km depth
- **Elements**: 9,214 triangular elements
- **Physics**: Rate-and-state friction with aging law
- **Duration**: 500 years simulation time

### Quick Start Guide

#### Step 1: Copy Example to Working Directory
```bash
# Create a working copy of example1
cp -r example1/ my_simulation/
cd my_simulation/
```

#### Step 2: Compile the Code
```bash
# Compile TriGreen for stiffness calculation
cd ../TriGreen/
./runcompile.sh

# Compile main simulation code
cd ../src/
./compile.sh
cd ../my_simulation/
```

#### Step 3: Calculate Stiffness Matrix
```bash
# Copy mesh file to TriGreen directory
cp triangular_mesh.gts ../TriGreen/

# Run stiffness calculation (single process for this example)
cd ../TriGreen/
mpirun -np 1 ./calc_trigreen

# Copy stiffness files back to example directory
cp trigreen_0.bin ../my_simulation/
cd ../my_simulation/
```

#### Step 4: Run the Simulation
```bash
# Run the earthquake cycle simulation
mpirun -np 1 ../src/3dtri_BP5 < parameter1.txt
```

### Example1 File Structure

```
example1/
├── parameter1.txt              # Main simulation parameters
├── triangular_mesh.gts         # Fault mesh geometry  
├── var-BP5_h1000.dat          # On-fault friction parameters
├── area-BP5_h1000.dat         # Element area data
├── profdp-BP5_h1000.dat       # Dip profile coordinates
├── profstrk-BP5_h1000.dat     # Strike profile coordinates
├── sub_stiff.sh               # SLURM script for stiffness calculation
└── sub_3dtri.sh               # SLURM script for main simulation
```

### Key Parameters in Example1

From `parameter1.txt`:
```bash
-BP5_h1000.dat              # File suffix for input files
result6/                    # Output directory  
5 9214 111 1 1 9 74 45     # Nab=5, Nt_all=9214, nprocs=1
31.50                       # Plate velocity: 31.5 mm/yr
500.0                       # Simulation time: 500 years
1.0 183.0 305.0            # Velocity thresholds (mm/s, mm/yr)
```

### Running on HPC Systems

#### Option 1: Interactive Mode
```bash
# Request compute node
salloc --nodes=1 --ntasks-per-node=1 --cpus-per-task=16 --time=2:00:00

# Set OpenMP threads
export OMP_NUM_THREADS=16
export OMP_PROC_BIND=close

# Run stiffness calculation
cd TriGreen/
mpirun -np 1 ./calc_trigreen

# Run simulation  
cd ../my_simulation/
mpirun -np 1 ../src/3dtri_BP5 < parameter1.txt
```

#### Option 2: Batch Submission
```bash
# Submit stiffness calculation
sbatch sub_stiff.sh

# Wait for completion, then submit main simulation
sbatch sub_3dtri.sh
```

### Expected Output Files

After successful completion, you should see:
```bash
result6/                        # Output directory
├── area-BP5_h1000.dat         # Updated area information
├── rupture-BP5_h1000.dat      # Rupture data
├── summary-BP5_h1000.dat      # Simulation summary
├── fltst_strk-*.dat           # Strike profiles
├── fltst_dip-*.dat            # Dip profiles  
├── timeseries_data_*.h5       # HDF5 time series (if enabled)
├── timeseries_data_*.xdmf     # Paraview visualization files
└── [various monitoring files]
```

### Performance Expectations

**Single Process (as configured):**
- **Stiffness calculation**: ~30-60 minutes
- **500-year simulation**: ~2-4 hours  
- **Memory usage**: ~8-12 GB
- **Disk space**: ~1-2 GB for outputs

**Scaling Options:**
```bash
# For faster execution, modify parameter1.txt:
# Change: 5 9214 111 1 1 9 74 45
# To:     5 9214 111 1 4 9 74 45  (4 processes)

# Then run with:
mpirun -np 4 ./calc_trigreen     # Stiffness
mpirun -np 4 ../src/3dtri_BP5    # Simulation
```

### Visualization

#### Option 1: Paraview (Recommended)
```bash
# Open XDMF files directly in Paraview
paraview result6/timeseries_data_-BP5_h1000.dat.xdmf
```

#### Option 2: MATLAB Analysis
```bash
# Use provided MATLAB scripts in Mesh/ directory
cd ../Mesh/
matlab -r "ReadInpFile; PlotVariable"
```

### Troubleshooting Example1

#### Common Issues:

**"TriGreen file not found":**
```bash
# Ensure stiffness files exist
ls trigreen_*.bin
# If missing, re-run stiffness calculation
```

**"Parameter file errors":**
```bash
# Check parameter1.txt format
# Ensure no extra spaces or missing lines
# Verify nprocs matches number of stiffness files
```

**Memory errors:**
```bash
# Reduce problem size or request more memory
# For SLURM: --mem-per-cpu=20G
```

**Slow performance:**
```bash
# Enable OpenMP
export OMP_NUM_THREADS=16
# Use multiple MPI processes
mpirun -np 4 ./3dtri_BP5
```

### Scientific Context

This example reproduces the SCEC SEAS BP5 benchmark, which models:
- **Long-term earthquake cycles** (~100-200 year recurrence)  
- **Interseismic loading** at plate velocity
- **Coseismic ruptures** with dynamic weakening
- **Postseismic slip** and stress relaxation

The results can be compared with other codes participating in the SCEC SEAS project for validation.

### Next Steps

After successfully running example1:
1. **Modify parameters** to explore different scenarios
2. **Scale up** to multiple processes for larger problems  
3. **Analyze results** using Paraview or MATLAB
4. **Create custom problems** using your own fault geometries

