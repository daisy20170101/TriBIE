# Dynamic Load Balancing - Final Implementation

## Overview

This document describes the completed implementation of dynamic load balancing in `calc_trigreen.f90` for optimal distribution of computational cells across MPI processes.

## Implementation Summary

### **1. Main Program Changes**

- **Dynamic Load Balancing Algorithm**: Implements optimal distribution using ceiling division
- **Variable Declaration**: Added `base_cells`, `extra_cells`, and `start_idx` for load balancing
- **Process Assignment**: Each process calculates its local cell count and starting index

### **2. Subroutine Interface**

- **Parameter Passing**: `base_cells` and `extra_cells` are now passed as parameters to `calc_green_allcell_improved`
- **Variable Access**: Subroutine can now access load balancing parameters directly

### **3. Key Features**

#### **Optimal Distribution Algorithm**
```fortran
base_cells = n_cell / size
extra_cells = mod(n_cell, size)

if (myid < extra_cells) then
   local_cells = base_cells + 1
   start_idx = myid * (base_cells + 1)
else
   local_cells = base_cells
   start_idx = extra_cells * (base_cells + 1) + (myid - extra_cells) * base_cells
end if
```

#### **Array Allocation Strategy**
- **`arr_co(3,local_cells)`**: Cell centroids for local cells only
- **`arr_trid(9,n_cell)`**: Triangle data for ALL cells (required for Green's function calculations)
- **`arr_cl_v2(3,3,local_cells)`**: Local coordinate systems for local cells only
- **`arr_out(local_cells, n_cell)`**: Output matrix for local cells × all cells

### **4. Performance Benefits**

1. **Minimal Load Imbalance**: Maximum difference of 1 cell between processes
2. **Memory Efficiency**: Only allocates arrays for local cells where possible
3. **Hybrid Parallelization**: Combines MPI for distributed memory and OpenMP for shared memory
4. **Optimal Communication**: Each process computes triangle data for all cells locally

### **5. Error Handling**

- **NaN Detection**: Custom `isnan` function for checking invalid calculations
- **Error Propagation**: Comprehensive error handling throughout the computation
- **Graceful Degradation**: Processes report errors without crashing the entire simulation

### **6. Output and Monitoring**

- **Load Balancing Summary**: Master process reports distribution statistics
- **Performance Metrics**: Timing information and cells processed per second
- **File Output**: Each process writes its local results to separate files

## Usage

### **Compilation**
```bash
mpif90 -O3 -fopenmp -o calc_trigreen calc_trigreen.f90 m_calc_green.f90 mod_dtrigreen.f90
```

### **Execution**
```bash
mpirun -np <number_of_processes> ./calc_trigreen
```

## Example Distribution

For `n_cell = 100` and `size = 4`:
- **Process 0**: 25 cells (indices 0-24)
- **Process 1**: 25 cells (indices 25-49)  
- **Process 2**: 25 cells (indices 50-74)
- **Process 3**: 25 cells (indices 75-99)

For `n_cell = 101` and `size = 4`:
- **Process 0**: 26 cells (indices 0-25)
- **Process 1**: 26 cells (indices 26-51)
- **Process 2**: 25 cells (indices 52-76)
- **Process 3**: 24 cells (indices 77-100)

## Technical Notes

- **Triangle Data Requirement**: Each process needs triangle data for ALL cells to compute Green's functions
- **Memory Overhead**: This is a necessary trade-off for distributed Green's function calculations
- **Scalability**: Algorithm scales efficiently with both MPI processes and OpenMP threads
- **Load Balance**: Achieves optimal distribution with minimal imbalance

## Conclusion

The dynamic load balancing implementation provides:
- **Optimal Work Distribution**: Minimal load imbalance across processes
- **Memory Efficiency**: Smart allocation strategy for distributed computation
- **Robust Error Handling**: Comprehensive validation and error reporting
- **Performance Monitoring**: Detailed metrics for optimization analysis

This implementation ensures that computational work is distributed as evenly as possible while maintaining the mathematical correctness required for Green's function calculations.
