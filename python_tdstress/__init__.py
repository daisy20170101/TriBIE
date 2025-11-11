"""
Python implementation of Triangular Dislocation stress/strain calculations.

Based on MATLAB code by Nikkhoo & Walter (2015).

Main functions:
--------------
tdstress_fs : Calculate stress/strain in elastic full-space
tdstress_hs : Calculate stress/strain in elastic half-space (complete implementation)

Reference:
---------
Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical,
artefact-free solution. Geophysical Journal International.
"""

from .tdstress_fs import tdstress_fs
from .tdstress_hs import tdstress_hs, tdstress_harfunc
from .td_utils import coord_trans, tens_trans, trimodefinder
from .ang_dislocation import ang_dis_strain, td_setup_s
from .ang_dislocation_fsc import ang_dis_strain_fsc
from .ang_setup_fsc import ang_setup_fsc_s

__version__ = "2.0.0"
__all__ = [
    'tdstress_fs',
    'tdstress_hs',
    'tdstress_harfunc',
    'coord_trans',
    'tens_trans',
    'trimodefinder',
    'ang_dis_strain',
    'td_setup_s',
    'ang_dis_strain_fsc',
    'ang_setup_fsc_s'
]
