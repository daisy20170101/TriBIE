# TriBIE_cycle
TriBIE_cycle is a Fortran90 parallel (MPI) code to simulation slow slip events and earthquake cycles in 3D curved fault geometry.

## Static Green's function (or stiffness) calculation
### Stiffness in homogeneous half-space
We use analytical solution for stress tensor in polygon elements embeded in homogeneous half-space.

### Using pylith for heterogeneous half-space
We use pylith to calculate stiffness for triangular elements in 3D curved slab geometry with either 1D or 3D velocity. 

## Compilation

``
$ mpiifort src/phy3d_module_non.f90 src/3dtri_cleanup.f90 -o 3dtri_sse
``

## Parameter files

## Simulation

