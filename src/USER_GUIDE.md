# TriBIE User Guide - Stiffness Computation and Cycling Simulation

## Overview
This guide explains how to use the TriBIE codebase for:
1. **Stiffness Matrix Computation** using `calc_trigreen.f90`
2. **Cycling Simulation** using `3dtri_BP5.f90` 
3. **Parameter Configuration** and file structures

## Workflow Overview

```
Mesh Generation → Stiffness Computation → Cycling Simulation → Analysis
     ↓                    ↓                    ↓              ↓
  triangular_mesh.gts  calc_trigreen.f90  3dtri_BP5.f90   Results
```

---

## 1. Stiffness Matrix Computation

### Prerequisites
- **Input Mesh**: `triangular_mesh.gts` file containing triangular fault elements
- **MPI Environment**: Multi-process execution environment
- **Compilation**: Use `TriGreen/runcompile.sh`

### Execution
```bash
cd TriGreen/
./runcompile.sh
mpirun -np <n_processes> ./calc_trigreen
```

### Output Files
- **`trigreen_<process_id>.bin`**: Stiffness matrix for each process
- **`position.bin`**: Element centroid positions
- **Console Output**: Distribution information and performance metrics

### Key Features
- **Dynamic Load Balancing**: Automatically distributes work unevenly across processes
- **SIMD Optimization**: Vectorized computation for improved performance
- **Memory Management**: Optimized allocation/deallocation strategies

---

## 2. Cycling Simulation

### Prerequisites
- **Stiffness Files**: `trigreen_<process_id>.bin` files from step 1
- **Parameter File**: `input/parameter1.txt` with simulation parameters
- **Compilation**: Use `src/compile.sh`

### Execution
```bash
cd src/
./compile.sh
cd ../input
mpirun -np <n_processes> ../src/3dtri_BP5
```

### Key Features
- **MPI_Scatterv**: Proper handling of uneven distributions
- **Rate-and-State Friction**: Physics-based fault behavior modeling
- **SIMD Optimization**: Vectorized physics calculations
- **Dynamic Load Balancing**: Compatible with calc_trigreen distribution

---

## 3. Parameter File Structure

### `parameter1.txt` Format
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

### Parameter Descriptions

#### **Array Dimensions**
- **`Nab`**: Number of temperature-dependent a-b profiles
- **`Nt_all`**: Total number of fault elements
- **`Nt`**: Local number of elements per process (auto-calculated)
- **`Lratio`**: Length ratio parameter
- **`nprocs`**: Number of MPI processes
- **`n_obv`**: Number of observation points
- **`np1`**: Number of strike profile points
- **`np2`**: Number of dip profile points

#### **Control Flags**
- **`Idin`**: Input mode flag
- **`Idout`**: Output mode flag
- **`Iprofile`**: Profile type selection
- **`Iperb`**: Perturbation flag
- **`Isnapshot`**: Snapshot output flag

#### **Physical Parameters**
- **`Vpl`**: Plate velocity (typically 1e-10 to 1e-9 m/s)
- **`tmax`**: Maximum simulation time in years

#### **Output Control**
- **`tint_out`**: Main output interval
- **`tmin_out`**: Minimum output interval
- **`tint_cos`**: Cosseismic output interval
- **`tint_sse`**: Slow slip event output interval

---

## 4. File Structure and Dependencies

### Required Input Files
```
src/
├── parameter1.txt           ! Main parameter file
├── triangular_mesh.gts     ! Fault mesh (if not using TriGreen)
└── [stiffness files]       ! From calc_trigreen or ssGreen format
```

### Generated Output Files
```
src/
├── area<jobname>           ! Area information
├── fltst_strk-*.dp+*      ! Strike-slip profiles
├── profstrk<jobname>       ! Strike profile data
├── profdp<jobname>         ! Dip profile data
└── [restart files]         ! If restart enabled
```

---

## 5. Example Workflow

### Step 1: Prepare Mesh
```bash
# Ensure triangular_mesh.gts exists in TriGreen/ directory
ls TriGreen/triangular_mesh.gts
```

### Step 2: Compute Stiffness
```bash
cd TriGreen/
./runcompile.sh
mpirun -np 8 ./calc_trigreen
# Generates: trigreen_0.bin, trigreen_1.bin, ..., trigreen_7.bin
```

### Step 3: Configure Parameters
```bash
cd src/
# Edit parameter1.txt with your simulation parameters
# Ensure nprocs matches the number of stiffness files
```

### Step 4: Run Simulation
```bash
./compile.sh
mpirun -np 8 ./3dtri_BP5
```

### Step 5: Analyze Results
```bash
# Check output files for results
ls -la area* fltst* prof*
```

---

## 6. Troubleshooting

### Common Issues

#### **MPI Communication Errors**
- **Symptom**: `MPI_ERR_TRUNCATE` or UCX crashes
- **Solution**: Ensure `nprocs` in parameter1.txt matches the number of stiffness files

#### **File Not Found Errors**
- **Symptom**: "TriGreen file not found" messages
- **Solution**: Verify stiffness files exist and are in the correct location

#### **Memory Issues**
- **Symptom**: Segmentation faults or allocation failures
- **Solution**: Reduce problem size or increase available memory

#### **Compilation Errors**
- **Symptom**: Compiler errors during build
- **Solution**: Ensure all dependencies are installed and compiler flags are correct

### Performance Optimization
- **SIMD**: Use `*_simd.f90` versions for better performance
- **Process Count**: Match MPI processes to available cores
- **Memory**: Monitor memory usage during large simulations

---

## 7. Advanced Features

### Dynamic Load Balancing
- Automatically handles uneven element distributions
- Optimizes work distribution across processes
- Compatible with irregular fault geometries

### SIMD Vectorization
- Automatic vectorization of computational loops
- Cache-aware memory access patterns
- Performance improvements on modern processors

### Restart Capability
- Save simulation state for later continuation
- Useful for long-running simulations
- Handles checkpoint/restart seamlessly

---

## 8. Output Analysis

### Key Output Files
1. **`area<jobname>`**: Element areas and properties
2. **`fltst_strk-*.dp+*`**: Strike-slip displacement profiles
3. **`profstrk<jobname>`**: Strike direction profiles
4. **`profdp<jobname>`**: Dip direction profiles

### Data Formats
- **Binary**: Most output files use unformatted binary format
- **Text**: Profile files use formatted text for easy parsing
- **Restart**: State files for simulation continuation

---

## 9. Best Practices

### Performance
- Use SIMD-optimized versions for production runs
- Match MPI processes to available hardware cores
- Monitor memory usage and adjust problem size accordingly

### Reliability
- Always verify input file formats and parameters
- Use restart capability for long simulations
- Check output files for expected results

### Maintenance
- Keep backup copies of working parameter sets
- Document any custom modifications
- Test with smaller problems before large-scale runs

---

## 10. Support and Resources

### Documentation
- This user guide
- Code comments and inline documentation
- Testing notes and examples

### Troubleshooting
- Check console output for error messages
- Verify file permissions and paths
- Ensure MPI environment is properly configured

### Performance Tuning
- Monitor timing output during execution
- Adjust block sizes for your specific hardware
- Profile memory usage patterns
