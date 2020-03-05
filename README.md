# SlowSlipSimulation
Model episodic slow slip events with 3D curved fault geometry


# Static Green's function (or stiffness) calculation

We could use pylith to calculate stiffness for triangular elements in 3D curved slab geometry with either 1D or 3D velocity. 
=======
## Compilation

### Method 1
``
$ scons
``

### Method 2

``
$ mpiifort src/phy3d_module_non.f90 src/3dtri_cleanup.f90 -o 3dtri_sse
``

