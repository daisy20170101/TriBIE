#! /usr/bin/python

# This is the configuration file for compiling Tri3DSSE

import os
import sys
#import libs


env=Environment(ENV={'PATH':os.environ['PATH']},LINK='mpiifort',F90='mpiifort')
#env.Tool('mpiifort')

sources=['3dtri_lock20_old.f90', 'phy3d_module_non.f90']

objs= env.Program(target='3dtri',source = ['src/phy3d_module_non.f90','src/3dtri_lock20_old.f90'])



#mpiifort phy3d_module_non.f90 3dtri_lock20_old.f90 -o 3dtri 
