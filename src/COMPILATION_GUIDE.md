# Compilation Guide for Improved Earthquake Simulation

## Overview

This guide explains how to compile the improved earthquake simulation code using either the provided bash script or Makefile.

## Prerequisites

### Required Software
- **MPI Implementation**: OpenMPI, MPICH, or Intel MPI
- **Fortran Compiler**: gfortran, ifort, or pgfortran (with MPI support)
- **Bash Shell**: For the compilation script (on Unix-like systems)

### Check Installation
```bash
# Check if mpif90 is available
which mpif90

# Check MPI version
mpif90 --version

# Check if you can run MPI programs
mpirun --version
```

## Compilation Methods

### Method 1: Using the Bash Script (Recommended)

The `compile.sh` script provides a comprehensive compilation process with error checking and colored output.

#### Basic Usage
```bash
# Make the script executable (first time only)
chmod +x compile.sh

# Compile the code
./compile.sh

# Compile and create symlink in parent directory
./compile.sh --symlink
```

#### What the Script Does
1. Checks for MPI Fortran compiler availability
2. Creates a `build/` directory
3. Compiles modules in dependency order
4. Compiles the main program
5. Links the executable
6. Provides compilation statistics
7. Optionally creates a symlink

#### Output
- **Executable**: `build/earthquake_simulation`
- **Object files**: `build/*.o`
- **Module files**: `build/*.mod`

### Method 2: Using the Makefile

The `Makefile` provides a traditional build system with dependency management.

#### Basic Usage
```bash
# Build the executable
make

# Clean build artifacts
make clean

# Rebuild from scratch
make rebuild

# Install (create symlink in parent directory)
make install

# Show available targets
make help
```

#### Makefile Targets
- `all` (default): Build the executable
- `clean`: Remove build artifacts
- `rebuild`: Clean and rebuild
- `install`: Create symlink in parent directory
- `uninstall`: Remove symlink
- `help`: Show help information

## Compilation Order

The modules are compiled in the following dependency order:

1. **`physical_constants.f90`** - No dependencies
2. **`simulation_parameters.f90`** - Depends on physical_constants
3. **`physical_variables.f90`** - Depends on physical_constants
4. **`io_parameters.f90`** - No dependencies
5. **`error_handling.f90`** - No dependencies
6. **`time_integration.f90`** - Depends on physical_constants
7. **`physics_equations.f90`** - Depends on physical_constants, physical_variables, simulation_parameters
8. **`mpi_utilities.f90`** - No dependencies
9. **`io_handling.f90`** - Depends on physical_constants, io_parameters
10. **`earthquake_simulation_improved.f90`** - Depends on all modules

## Compiler Flags

### Default Flags
```bash
-O2          # Optimization level 2
-g           # Include debug information
-Wall        # Enable all warnings
-Wextra      # Enable extra warnings
-std=f2008   # Use Fortran 2008 standard
```

### Customizing Flags
You can modify the compiler flags in either:

#### Bash Script
Edit the `COMPILER_FLAGS` variable in `compile.sh`:
```bash
COMPILER_FLAGS="-O3 -g -Wall -Wextra -std=f2008 -fopenmp"
```

#### Makefile
Edit the `FFLAGS` variable in `Makefile`:
```makefile
FFLAGS = -O3 -g -Wall -Wextra -std=f2008 -fopenmp
```

## Troubleshooting

### Common Issues

#### 1. MPI Compiler Not Found
```bash
Error: MPI Fortran compiler 'mpif90' not found!
```
**Solution**: Install an MPI implementation or ensure it's in your PATH.

#### 2. Compilation Errors
Check the error messages for:
- Missing dependencies
- Syntax errors
- Module import issues

#### 3. Linking Errors
Ensure all object files are created and the linking order is correct.

### Debug Mode
For debugging, use more verbose flags:
```bash
# In compile.sh or Makefile
FFLAGS="-O0 -g -Wall -Wextra -std=f2008 -fcheck=all -fbacktrace"
```

### Verbose Output
The bash script provides detailed output. For Makefile, use:
```bash
make VERBOSE=1
```

## Running the Simulation

After successful compilation:

```bash
# Navigate to build directory
cd build

# Run with MPI (example: 4 processes)
mpirun -np 4 ./earthquake_simulation

# Or if you used --symlink or make install
mpirun -np 4 ../earthquake_simulation
```

## Performance Optimization

### Compiler Optimizations
- **`-O2`**: Good balance of optimization and compilation speed
- **`-O3`**: Maximum optimization (may increase compilation time)
- **`-march=native`**: Optimize for current CPU architecture

### MPI Optimizations
- **`-O2`**: Standard optimization level
- **`-O3`**: Aggressive optimization
- **`-ffast-math`**: Fast math operations (may affect precision)

### Example Optimized Flags
```bash
FFLAGS="-O3 -g -march=native -ffast-math -funroll-loops"
```

## Platform-Specific Notes

### Linux
- Works with most MPI implementations
- gfortran is commonly available

### macOS
- OpenMPI works well
- May need to install gfortran separately

### Windows
- Use WSL (Windows Subsystem for Linux)
- Or use Intel MPI with Intel Fortran

## Support

If you encounter compilation issues:

1. Check the error messages carefully
2. Verify MPI installation
3. Ensure all source files are present
4. Check compiler compatibility
5. Review the dependency order

The compilation scripts include error checking and will provide helpful information about what went wrong.
