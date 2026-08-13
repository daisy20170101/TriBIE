# Testing Notes - MPI_Scatterv Implementation

## Test Date
December 19, 2024

## Test Environment
- **Model**: SEAS BP3 fault model
- **Implementation**: MPI_Scatterv for uneven distribution handling
- **Files Modified**: 
  - `3dtri_BP5.f90` - Original file with MPI_Scatterv implementation
  - `3dtri_BP5_simd.f90` - SIMD-optimized version with MPI_Scatterv implementation

## Test Results
✅ **SUCCESSFUL** - The new MPI_Scatterv implementation has been tested with the SEAS BP3 fault model and is working correctly.

## What Was Tested
1. **Uneven Distribution Handling**: Verified that processes with different `local_cells` values can properly receive data
2. **MPI Communication**: Confirmed no more `MPI_ERR_TRUNCATE` errors
3. **Array Sizing**: Validated that arrays are correctly allocated with `local_cells` instead of `Nt`
4. **Memory Management**: Checked proper allocation/deallocation of `sendcounts` and `displs` arrays

## Technical Details
- **Before**: Used `MPI_Scatter` which requires uniform distribution (caused truncation errors)
- **After**: Implemented `MPI_Scatterv` which handles uneven distributions properly
- **Array Allocation**: Fixed from `Nt` to `local_cells` to match actual data size per process
- **Distribution Algorithm**: Compatible with `calc_trigreen.f90` dynamic load balancing

## Status
**READY FOR PRODUCTION USE** - The MPI_Scatterv implementation resolves the UCX communication crashes and message truncation errors that were occurring with the SEAS BP3 fault model.
