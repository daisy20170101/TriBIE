# Dynamic Load Balancing in TriGreen

## Overview

This implementation provides optimal distribution of computational cells across MPI processes when the total number of cells (`n_cell`) cannot be evenly divided by the number of processes (`size`).

## How It Works

### 1. **Optimal Distribution Algorithm**
- **Base cells per process**: `base_cells = n_cell / size`
- **Extra cells to distribute**: `extra_cells = mod(n_cell, size)`
- **Distribution strategy**: First `extra_cells` processes get `base_cells + 1` cells, remaining processes get `base_cells` cells

### 2. **Process Assignment**
```fortran
if (myid < extra_cells) then
   Nt = base_cells + 1
   start_idx = myid * (base_cells + 1)
else
   Nt = base_cells
   start_idx = extra_cells * (base_cells + 1) + (myid - extra_cells) * base_cells
end if
```

### 3. **Communication Pattern**
- **Master process (Rank 0)**: Calculates distribution and sends individual messages to each worker
- **Worker processes**: Receive their specific cell count and starting index
- **Module variables**: Share distribution parameters between main program and subroutines

## Example Distributions

### Case 1: Even Distribution
- `n_cell = 100`, `size = 4`
- Each process gets 25 cells
- No load imbalance

### Case 2: Uneven Distribution
- `n_cell = 101`, `size = 4`
- **Process 0**: 26 cells (indices 0-25)
- **Process 1**: 26 cells (indices 26-51)
- **Process 2**: 25 cells (indices 52-76)
- **Process 3**: 25 cells (indices 77-101)
- **Load imbalance**: 1 cell (excellent balance!)

### Case 3: Large Uneven Distribution
- `n_cell = 103`, `size = 4`
- **Process 0**: 26 cells (indices 0-25)
- **Process 1**: 26 cells (indices 26-51)
- **Process 2**: 26 cells (indices 52-77)
- **Process 3**: 25 cells (indices 78-103)
- **Load imbalance**: 1 cell

## Benefits

1. **No Lost Cells**: All cells are processed
2. **Minimal Imbalance**: Maximum difference of 1 cell between processes
3. **Efficient Communication**: Individual messages instead of broadcasts
4. **Scalable**: Works for any number of cells and processes
5. **Transparent**: Clear reporting of distribution and performance

## Performance Monitoring

The implementation includes:
- Distribution summary at startup
- Process-specific cell count reporting
- Performance timing and statistics
- Load balancing efficiency metrics

## Usage

1. **Compile**: Use the provided compilation script
2. **Run**: Execute with desired number of MPI processes
3. **Monitor**: Check output for distribution information and performance metrics

## Files Modified

- `calc_trigreen.f90`: Main program with dynamic load balancing
- `m_calc_green.f90`: Module with shared load balancing variables
- `test_load_balancing.f90`: Test program demonstrating the algorithm
- `compile_test.sh`: Script to compile and test the load balancing

## Future Enhancements

- Adaptive load balancing based on process performance
- Dynamic redistribution during computation
- Load balancing for heterogeneous systems
- Integration with OpenMP for hybrid parallelization
