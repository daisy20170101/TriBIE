# SlowSlipSimulation
Model episodic slow slip events with 3D curved fault geometry

## Compilation

### Method 1
``
$ scons
``

### Method 2

``
$ mpiifort src/phy3d_module_non.f90 src/3dtri_cleanup.f90 -o 3dtri_sse
``

