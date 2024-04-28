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

file list:
parameters.txt: key initial configuration of simulation
on-fault parameters: var-BP5_h500_140_60.dat
observation files: profdp-BP5_h500_140_60.dat, profstrk-BP5_h500_140_60.dat
mesh element area: area-BP5_h500_140_60.dat
mesh file: fault_h500_140_60.gts


## Simulation

![varBP5](https://github.com/daisy20170101/TriBIE/assets/33549997/c0b43d1b-777a-48e0-bda4-72c7a9b0e95e)
Figure: Mapview of on-fault distribution of fault key parameters in BP5 example.

## Results
![bp5_slip_h1000_140_60](https://github.com/daisy20170101/TriBIE/assets/33549997/b5e8804c-297d-4bbd-b8cf-c75a64f6bea9)
Figure: Cumulative slip along strike (left) and along downdip (right) in the first 800 modeling years. Coseismic slip in red while interseismic in blue. 



