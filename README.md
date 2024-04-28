# TriBIE
TriBIE is a Fortran90 parallel (MPI) code to simulation slow slip events and earthquake cycles in an arbitrary curved fault buried in the half-space medium. TriBIE has been verified in the SCEC Sequences of Earthquakes and Aseismic Slip Project (https://strike.scec.org/cvws/seas/).

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

![varBP5](https://github.com/daisy20170101/TriBIE/assets/33549997/c0b43d1b-777a-48e0-bda4-72c7a9b0e95e)
Figure: Mapview of on-fault distribution of fault key parameters in BP5 example.

