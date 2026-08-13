# Direct TriGreen Integration for 3dtri_BP5.f90

## Overview

Since the `trigreen_<process_id>.bin` files from `calc_trigreen.f90` are already in the same format as the `ssGreen_<process_id>.bin` files expected by `3dtri_BP5.f90`, we can integrate them directly by modifying the loading mechanism in `3dtri_BP5.f90`.

## Key Insight

**No file conversion is needed!** The TriGreen output files are already compatible with the stiffness matrix format. We just need to:

1. **Calculate the correct `local_cells`** for each process using the same dynamic load balancing as `calc_trigreen.f90`
2. **Modify the file loading** to use `trigreen_<process_id>.bin` instead of `ssGreen_<process_id>.bin`
3. **Ensure the array allocation** matches the actual number of local cells

## Required Changes

### 1. Add Dynamic Load Balancing Variables

Add these variables after the MPI definitions:

```fortran
! Dynamic load balancing variables (compatible with calc_trigreen.f90)
integer :: base_cells, extra_cells, local_cells, start_idx
logical :: use_trigreen_format = .true.  ! Set to .true. to use TriGreen files
```

### 2. Implement Dynamic Load Balancing

Replace the existing load balancing logic with:

```fortran
! MODIFICATION: Implement dynamic load balancing compatible with calc_trigreen.f90
if (use_trigreen_format) then
   ! Calculate optimal distribution using the same algorithm as calc_trigreen.f90
   base_cells = Nt_all / size
   extra_cells = mod(Nt_all, size)
   
   if (myid < extra_cells) then
      local_cells = base_cells + 1
      start_idx = myid * (base_cells + 1)
   else
      local_cells = base_cells
      start_idx = extra_cells * (base_cells + 1) + (myid - extra_cells) * base_cells
   end if
   
   ! Override Nt with the actual local cells for this process
   Nt = local_cells
   
   if (myid == master) then
      write(*,*) 'Using TriGreen format with dynamic load balancing'
      write(*,*) 'Base cells per process:', base_cells
      write(*,*) 'Extra cells distributed:', extra_cells
      write(*,*) 'Total processes:', size
   end if
   
   write(*,*) 'Process', myid, 'gets', local_cells, 'cells starting from index', start_idx
else
   ! Original logic for even distribution
   if(mod(Nt_all,nprocs)/=0)then
      write(*,*)'Nd_all must be integer*nprocs. Change nprocs!'
      STOP
   else
      write(*,*)'Each cpu calculates',Nt_all/nprocs,'cells'
   end if
end if
```

### 3. Update Array Allocation

Change the stiffness matrix allocation from:

```fortran
ALLOCATE (stiff(Nt,Nt_all),stiff2(Nt,Nt_all))
```

To:

```fortran
ALLOCATE (stiff(local_cells,Nt_all),stiff2(local_cells,Nt_all))
```

### 4. Modify File Loading

Change the file opening from:

```fortran
open(5, file=trim(stiffname)//'ssGreen_'//trim(adjustl(cTemp))//'.bin',form='unformatted',access='stream',buffered='yes')
```

To:

```fortran
if (use_trigreen_format) then
   ! Load TriGreen format files
   open(5, file=trim(stiffname)//'trigreen_'//trim(adjustl(cTemp))//'.bin',form='unformatted',access='stream',buffered='yes')
   write(*,*) 'Loading TriGreen file: trigreen_', trim(adjustl(cTemp)), '.bin'
else
   ! Load original ssGreen format files
   open(5, file=trim(stiffname)//'ssGreen_'//trim(adjustl(cTemp))//'.bin',form='unformatted',access='stream',buffered='yes')
   write(*,*) 'Loading ssGreen file: ssGreen_', trim(adjustl(cTemp)), '.bin'
end if
```

### 5. Update Loops and MPI Operations

Change the stiffness matrix reading loop from:

```fortran
do i=1,Nt !! observe
```

To:

```fortran
do i=1,local_cells !! observe (now using local_cells instead of Nt)
```

Update MPI_Scatter calls from:

```fortran
call MPI_Scatter(cca_all,Nt,MPI_Real8,cca,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
```

To:

```fortran
call MPI_Scatter(cca_all,local_cells,MPI_Real8,cca,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
```

## Complete Workflow

### Step 1: Run TriGreen
```bash
cd TriGreen
mpirun -np 4 ./calc_trigreen
```

This creates:
- `trigreen_0.bin`, `trigreen_1.bin`, `trigreen_2.bin`, `trigreen_3.bin`
- `position.bin`

### Step 2: Apply the Patch to 3dtri_BP5.f90
Use the provided patch file `3dtri_BP5_trigreen_patch.txt` to modify your `3dtri_BP5.f90`.

### Step 3: Run 3dtri_BP5.f90
```bash
cd ../src
mpirun -np 4 ./3dtri_BP5
```

## Benefits of Direct Integration

1. **No file conversion needed** - Direct use of TriGreen output
2. **Consistent load balancing** - Same algorithm used in both programs
3. **Memory efficient** - Arrays sized correctly for actual local cells
4. **Backward compatible** - Can switch between formats using the flag
5. **Performance optimized** - Maintains all existing optimizations

## Configuration

### Enable TriGreen Format
Set the flag in the code:
```fortran
logical :: use_trigreen_format = .true.  ! Set to .true. to use TriGreen files
```

### Disable TriGreen Format (Use Original)
Set the flag to:
```fortran
logical :: use_trigreen_format = .false.  ! Use original ssGreen format
```

## File Structure

With TriGreen format enabled, the program expects:
```
<stiffname>/
├── trigreen_0.bin          # Green's functions for process 0
├── trigreen_1.bin          # Green's functions for process 1
├── trigreen_2.bin          # Green's functions for process 2
├── trigreen_3.bin          # Green's functions for process 3
├── position.bin             # Cell positions
└── (other required files)
```

## Troubleshooting

### Common Issues

1. **File not found errors**
   - Ensure TriGreen output files exist in the specified directory
   - Check that `stiffname` in `parameter1.txt` points to the correct directory

2. **Array bounds errors**
   - Verify that `local_cells` is calculated correctly
   - Check that all arrays are allocated with the correct dimensions

3. **MPI scatter errors**
   - Ensure `local_cells` is used consistently in all MPI operations
   - Verify that the total number of cells matches across all processes

### Verification Steps

1. **Check load balancing output** - Each process should report its assigned cells
2. **Verify file loading** - Each process should report loading its TriGreen file
3. **Check array dimensions** - Stiffness matrix should be `local_cells × Nt_all`

## Summary

This direct integration approach eliminates the need for file conversion while maintaining full compatibility with the existing `3dtri_BP5.f90` codebase. The key is implementing the same dynamic load balancing algorithm used in `calc_trigreen.f90` to ensure consistent cell distribution across processes.
