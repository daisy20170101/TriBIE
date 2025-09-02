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
<Nab> <Nt_all> <Nt> <Lratio> <nprocs> <n_obv> <np1> <np2>  ! Array dimensions
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

