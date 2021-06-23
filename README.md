# TriBIE
TriBIE is a Fortran90 parallel (MPI) code to simulation slow slip events and earthquake cycles in an arbitrary curved fault buried in the half-space medium. TriBIE has been verified in the SCEC Sequences of Earthquakes and Aseismic Slip Project (https://strike.scec.org/cvws/seas/).

## Static Green's function (or stiffness) calculation
### Stiffness in homogeneous half-space
We use analytical solution for stress tensor in polygon elements embeded in homogeneous half-space.

### Using pylith for heterogeneous half-space
We use pylith to calculate stiffness for triangular elements in 3D curved slab geometry with either 1D or 3D velocity. 

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

## Simulation

