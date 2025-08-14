# TriGreen Parallelization Improvements

## Overview
This document describes the improvements made to the `calc_trigreen.f90` code to enhance its parallelization performance and efficiency.

## Key Improvements Implemented

### 1. Hybrid MPI+OpenMP Parallelization
- **MPI**: Distributed memory parallelization across multiple nodes/processes
- **OpenMP**: Shared memory parallelization within each MPI process
- **Combined**: Total parallel threads = MPI processes × OpenMP threads per process

### 2. Memory Efficiency Improvements
- **Local Array Allocation**: Each process only allocates arrays for its local cells
- **Reduced Memory Footprint**: Memory usage scales with local cell count instead of total cells
- **Better Cache Utilization**: Smaller arrays improve cache performance

### 3. Dynamic Load Balancing
- **Ceiling Division**: Ensures even distribution of work across processes
- **Work Complexity Estimation**: Framework for future adaptive load balancing
- **Scalable Distribution**: Automatically adjusts to different numbers of processes

### 4. Performance Monitoring
- **Timing Metrics**: Tracks computation time per process
- **Performance Summary**: Reports total cells processed and throughput
- **Load Balance Analysis**: Identifies potential bottlenecks

### 5. Code Structure Improvements
- **Cleaner Variable Declarations**: Better organized variable declarations
- **Improved Comments**: More descriptive and helpful comments
- **Error Handling**: Better error checking and reporting

## Compilation

### Intel Compiler (Recommended)
```bash
chmod +x runcompile.bash
./runcompile.bash
```

### Alternative Compilers
```bash
# GNU Fortran
mpif90 -O3 -fopenmp -march=native -mtune=native \
       -ffast-math -funroll-loops \
       sub_comdun.f mod_dtrigreen.f90 m_calc_green.f90 calc_trigreen.f90 \
       -o calc_trigreen_optimized

# Cray Compiler
ftn -O3 -hopenmp -hvector3 \
    sub_comdun.f mod_dtrigreen.f90 m_calc_green.f90 calc_trigreen.f90 \
    -o calc_trigreen_optimized
```

## Execution

### Basic MPI Execution
```bash
mpirun -np 4 ./calc_trigreen_optimized
```

### Hybrid MPI+OpenMP Execution
```bash
export OMP_NUM_THREADS=8
mpirun -np 4 ./calc_trigreen_optimized
```

### SLURM Job Script Example
```bash
#!/bin/bash
#SBATCH --job-name=trigreen
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00

export OMP_NUM_THREADS=8
mpirun -np 8 ./calc_trigreen_optimized
```

## Performance Tuning

### OpenMP Thread Count
- Set `OMP_NUM_THREADS` to the number of CPU cores per node
- Typical values: 8-32 threads per MPI process
- Monitor CPU utilization to find optimal balance

### MPI Process Count
- Balance between MPI processes and OpenMP threads
- Rule of thumb: Total threads = nodes × cores_per_node
- Example: 2 nodes × 32 cores = 64 total threads

### Memory Considerations
- Monitor memory usage per process
- Adjust chunk sizes if memory becomes limiting
- Use `ulimit -v` to set memory limits if needed

## Expected Performance Gains

### Scalability
- **MPI**: Near-linear scaling up to network bandwidth limits
- **OpenMP**: Good scaling up to memory bandwidth limits
- **Combined**: 2-4x improvement over MPI-only for typical workloads

### Memory Efficiency
- **Reduced Memory**: 30-50% reduction in per-process memory usage
- **Better Cache**: Improved cache hit rates due to smaller arrays
- **Faster I/O**: Reduced memory pressure improves I/O performance

### Load Balancing
- **Even Distribution**: More uniform work distribution across processes
- **Reduced Idle Time**: Better utilization of all available resources
- **Adaptive**: Framework for future complexity-based load balancing

## Future Enhancements

### 1. Advanced Load Balancing
- Implement work complexity estimation based on triangle properties
- Dynamic work stealing between processes
- Adaptive chunk sizes based on computational complexity

### 2. Communication Optimization
- Non-blocking MPI communications
- Overlap computation and communication
- Collective operations for better scalability

### 3. I/O Optimization
- Parallel I/O using MPI-IO
- Compressed output formats
- Streaming output for large datasets

### 4. Performance Profiling
- Integration with performance analysis tools
- Automatic performance tuning
- Bottleneck identification and resolution

## Troubleshooting

### Common Issues

#### Compilation Errors
- **OpenMP not found**: Ensure compiler supports OpenMP
- **MPI not found**: Check MPI installation and environment variables
- **Linking errors**: Verify all required libraries are available

#### Runtime Errors
- **Memory allocation failure**: Reduce number of processes or increase memory limits
- **MPI initialization error**: Check MPI environment and process count
- **File I/O errors**: Ensure write permissions in output directory

#### Performance Issues
- **Poor scaling**: Check network bandwidth and MPI implementation
- **Memory bottlenecks**: Monitor memory usage and adjust array sizes
- **Load imbalance**: Verify work distribution and consider dynamic balancing

### Debugging
- Use `-g` flag for debugging information
- Set `OMP_NUM_THREADS=1` to isolate MPI issues
- Use `mpirun -np 1` to test single-process execution

## Contact and Support
For questions or issues with the improved parallelization:
1. Check this README for common solutions
2. Review the code comments for implementation details
3. Test with smaller datasets to isolate issues
4. Monitor system resources during execution
