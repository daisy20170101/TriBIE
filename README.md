# TriBIE
TriBIE is a parallelized Fortran90 (OpenMPI) code for simulating episodic slow slip events (SSEs) and earthquake cycles in an arbitrary curved fault buried in half-space medium. The quasi-dynamic formulation of shear traction and fault slip [Rice 1993] is implemented in spatial domain and the 5th-order Runge-Kutta method with adaptive time step integration is used in time domain. TriBIE has been used in many geophysical applications recently. TriBIE has also been verified in the SCEC Sequences of Earthquakes and Aseismic Slip Project (https://strike.scec.org/cvws/seas/). 

TriBIE is developed and currently maintaned by [Duo Li](https://daisy20170101.github.io/duo_li.github.io/), [Yajing Liu](https://www.mcgill.ca/eps/liu), [Hongyu Yu](https://scholar.google.com/citations?user=XzF1Ed0AAAAJ&hl=en), [Lei Zhang](https://www.eq-igl.ac.cn/fgzc/info/2019/20831.html) and [Andrea Perez-Silva](https://www.researchgate.net/scientific-contributions/Andrea-Perez-Silva-2190753609). Please feel free to contact us for anything about the code!

![image](https://github.com/daisy20170101/TriBIE/blob/master/post_precessing/SSE_CAS_web.png)
Figure 1. Numerical model of the curved subduction fault in Cascadia (left) and the simulated SSE slip in two episodes [Li and Liu, 2016; 2017].

## Static Green's function (or stiffness) calculation
### Stiffness in homogeneous half-space medium
We use analytical solution for stress tensor in polygon elements embeded in homogeneous half-space based on Comninou and Dundurs [1975] and Stuart [1997].
We fix the singularities using the method presented in Nikkhoo and Walter [2015, GJI]. The triangular mesh of fault interface can be created by [Trelis or Cubit](https://coreform.com/products/coreform-cubit/). 

![image](https://github.com/daisy20170101/TriBIE/blob/master/post_precessing/GuerreroSlabMesh.png)
Figure 2. Fault geometry and mesh created for SSE model in Guerrero, Mexico.

Please contact us if you need scripts for creating Trelis meshing file. 
    
### Using pylith for heterogeneous half-space medium
We can use [Pylith](https://geodynamics.org/resources/pylith) to calculate stiffness for triangular elements in 3D curved slab geometry with either 1D or 3D velocity. 

![image](https://github.com/daisy20170101/TriBIE/blob/master/post_precessing/CAS_3d2.png)
Figure 3. Mesh and 3D velocity model of Cascadia subduction used in PyLith simulation. In Pylith, only traction changes on the slab fault is needed.

Please contact us if you need help with calculating stiffness using PyLith. 

## Compilation

``
$ mpiifort src/phy3d_module_non.f90 src/3dtriBIE_v1.f90 -o 3dtriBIE
``

## Parameter files

Line 1: jobname (appendix)

Line 2: foldername

Line 3: stiffname

Line 4: restartname

Line 5: profile

Line 6: Number of points in a-b profile, total element,element/process,num. of procs,h*

Line 7: Idin,Idout,Iprofile,Iperb,Isnapshot

Line 8: factor1,factor2,factor3,factor4

Line 9: Vpl in mm/yr

Line 10: up-limit of locking depth, bottom-limit of locking depth

Line 11: \bar{sigma} in locking depth, \bar{sigma} in SSE,d_c in SSE depth

Line 12: simulation time in year

Line 13: tslip_ave,tslip_aveint

Line 14: tint_out,tmin_out,tint_cos

Line 15: vcos,vsse1,vsse2

Line 16: tssestart,tsseend

Line 17: nmv,nas,ncos,nsse

Line 18: s1,s2,s3,s4,s5,s6,s7,s8


Example:

    -test.dat
    ./out/
    ../Code1/TriGreen/
    ./out/out0-test.dat
    ../Code1/profiles/GuerreroT.dat
    5 90496 1414 64 122      !Nab,Nt_all,Nt,nprocs,hnucl
    1 1 2 1 0                !Idin, Idout Iprofile, Iperb,Isnapshot
    2.0  0.98  1.03 1.0      ! (factor1 factor2 factor3 factor4)
    61.0                     !    Vpl (mm/yr)
    20.0 45.0                !xilock1 xilock2 (SSE depth)
    500.0 25.0 10.0862       !sigmadiff, sefffix (bar), Lffix (mm)
    200.0                    ! tmax
    0.0 5.0                  !tslip_ave  tslip_aveint
    100.0 1.0 1.585d-8       !tint_out tmin_out tint_cos (3.17d-8 is 1sec)
    5.0 61.0 305.0           !vcos (mm/s), vsse1(3Vpl) vsse2 (5Vpl, mm/yr)
    0.0 2000.0                                           ! sse start and end times to record moment
    500 30 30 30                                         ! nmv,nas,ncos,nsse
    72081 15189 5372 4900 3138 65037 18071 20293         ! output element IDs (s1-s8)


Please contact us or post an issue if anything unclear!

## Publications

Jiang, J., B. Erickson, V. Lambert, J.-P. Ampuero, R. Ando, S. Barbot, C. Cattania, L. D. Zilio, B. Duan, E. M. D. . and et al. "Community-Driven Code Comparisons for Three-Dimensional Dynamic Modeling of Sequences of Earthquakes and Aseismic Slip (SEAS)." JGR: Solid Earth

Perez-Silva, A., D. Li, A.-A. Gabriel and Y. Kaneko (2021). "3D Modeling of Long-Term Slow Slip Events Along the Flat-Slab Segment in the Guerrero Seismic Gap, Mexico." Geophysical Research Letters 48(13).

Perez-Silva, A., Y. Kaneko, M. Savage, L. Wallace, D. Li and C. Williams (2022). "Segmentation of Shallow Slow Slip Events at the Hikurangi Subduction Zone Explained by Along-Strike Changes in Fault Geometry and Plate Convergence Rates." Journal of Geophysical Research: Solid Earth 127(1): e2021JB022913.

Li, H., M. Wei, D. Li, Y. Liu, Y. Kim and S. Zhou "Segmentation of Slow Slip Events in South Central Alaska Possibly Controlled by a Subducted Oceanic Plateau." Journal of Geophysical Research: Solid Earth.

H Yu, Y Liu, H Yang, J Ning , Modeling earthquake sequences along the Manila subduction zone: Effects of three-dimensional fault geometry. Tectonophysics, 2018

Li, D. and Y. Liu (2016). "Spatiotemporal evolution of slow slip events in a nonplanar fault model for northern Cascadia subduction zone." Journal of Geophysical Research: Solid Earth 121(9): 6828-6845.

Li, D. and Y. Liu (2017). "Modeling slow-slip segmentation in Cascadia subduction zone constrained by tremor locations and gravity anomalies." Journal of Geophysical Research: Solid Earth 122: 3138–3157.

