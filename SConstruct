#! /usr/bin/python

# This is the configuration file for compiling TriBIE

import os

mpi_fc = os.environ.get('MPI_FC', os.environ.get('MPIFC', 'mpifort'))

env = Environment(ENV={'PATH': os.environ['PATH']}, LINK=mpi_fc, F90=mpi_fc, FORTRAN=mpi_fc)

sources = ['src/phy3d_module_non.f90', 'src/3dtriBIE_v1.f90']

program = env.Program(target='3dtriBIE', source=sources)



#mpiifort phy3d_module_non.f90 3dtri_lock20_old.f90 -o 3dtri 
