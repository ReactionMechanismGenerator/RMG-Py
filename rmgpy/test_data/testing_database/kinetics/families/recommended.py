#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This file contains multiple sets of suggested kinetics families for various
systems of interest. They can be used by including the name of a set in the
kineticsFamilies part of the input file. Multiple sets can be specified at the
same time, and union of them will be loaded. These sets can also be specified
along with individual families. Custom sets can be easily defined in this file
and immediately used in input files without any additional changes.
"""

default = {
    'Disproportionation',
    'H_Abstraction',
    'R_Addition_MultipleBond',
    'R_Recombination',
}

pah = {
    '1,2_shiftC',
    '6_membered_central_C-C_shift',
    'Intra_R_Add_Exo_scission',
    'Intra_ene_reaction',
}
