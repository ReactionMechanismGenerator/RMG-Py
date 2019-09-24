#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
This script generates a video showing the flux diagram for a given reaction
model as it evolves in time. It takes as its arguments the path to an RMG-Py
input file corresponding to a job that has already been run and the
corresponding Chemkin mechanism and RMG species dictionary files. If a folder
of species images is available, it can be passed as an optional argument. A
Chemkin output file can also be passed as an optional positional argument.
"""

import argparse
import os

from rmgpy.tools.fluxdiagram import create_flux_diagram


################################################################################

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('input', metavar='INPUT', type=str, help='RMG input file')
    parser.add_argument('chemkin', metavar='CHEMKIN', type=str, help='Chemkin file')
    parser.add_argument('dictionary', metavar='DICTIONARY', type=str, help='RMG dictionary file')
    parser.add_argument('chemkinOutput', metavar='CHEMKIN_OUTPUT', type=str, nargs='?', default=None,
                        help='Chemkin output file')
    parser.add_argument('--java', action='store_true', help='process RMG-Java model')
    parser.add_argument('--no-dlim', dest='dlim', action='store_false', help='Turn off diffusion-limited rates')
    parser.add_argument('-s', '--species', metavar='DIR', type=str, help='Path to folder containing species images')
    parser.add_argument('-f', '--foreign', dest='checkDuplicates', action='store_true',
                        help='Not an RMG generated Chemkin file (will be checked for duplicates)')
    parser.add_argument('-n', '--maxnode', metavar='N', type=int, help='Maximum number of nodes to show in diagram')
    parser.add_argument('-e', '--maxedge', metavar='N', type=int, help='Maximum number of edges to show in diagram')
    parser.add_argument('-c', '--conctol', metavar='TOL', type=float, help='Lowest fractional concentration to show')
    parser.add_argument('-r', '--ratetol', metavar='TOL', type=float, help='Lowest fractional species rate to show')
    parser.add_argument('-t', '--tstep', metavar='S', type=float,
                        help='Multiplicative factor to use between consecutive time points')
    parser.add_argument('--centralSpecies', metavar='s1,s2,...', type=lambda s: [int(idx) for idx in s.split(',')],
                        help='List of indices of central species')
    parser.add_argument('--rad', metavar='R', type=int, help='Graph radius around a central species (only useful if'
                                                             ' not using super)')
    parser.add_argument('--centralReactionCount', metavar='N', type=int, default=1,
                        help='Maximum number of reactions to show from each central species (default = 1).'
                             ' If rad > 1, then this is the number of reactions from every species')
    parser.add_argument('--super', action='store_true', help='Superimpose central species onto normal flux diagram to'
                                                             ' ensure that they appear in diagram (might result in more'
                                                             ' nodes and edges than given by maxnode and maxedge)')
    parser.add_argument('--saveStates', action='store_true', help='Save simulation states to disk')
    parser.add_argument('--readStates', action='store_true', help='Read simulation states from disk')

    args = parser.parse_args()

    input_file = os.path.abspath(args.input)
    chemkin_file = os.path.abspath(args.chemkin)
    dict_file = os.path.abspath(args.dictionary)
    species_path = os.path.abspath(args.species) if args.species is not None else None
    chemkin_output = os.path.abspath(args.chemkinOutput) if args.chemkinOutput is not None else ''
    use_java = args.java
    dflag = args.dlim
    check_duplicates = args.checkDuplicates
    central_species_list = args.centralSpecies
    superimpose = args.super
    save_states = args.saveStates
    read_states = args.readStates

    keys = ('max_node_count',
            'max_edge_count',
            'concentration_tol',
            'species_rate_tol',
            'radius',
            'central_reaction_count',
            'time_step')
    vals = (args.maxnode, args.maxedge, args.conctol, args.ratetol, args.rad, args.central_reaction_count, args.tstep)
    settings = {k: v for k, v in zip(keys, vals) if v is not None}

    return (input_file,
            chemkin_file,
            dict_file,
            species_path,
            chemkin_output,
            use_java,
            dflag,
            check_duplicates,
            settings,
            central_species_list,
            superimpose,
            save_states,
            read_states)


def main():
    (input_file,
     chemkin_file,
     dict_file,
     species_path,
     chemkin_output,
     use_java,
     dflag,
     check_duplicates,
     settings,
     central_species_list,
     superimpose,
     save_states,
     read_states) = parse_arguments()

    create_flux_diagram(input_file, chemkin_file, dict_file, species_path=species_path, java=use_java,
                        settings=settings, chemkin_output=chemkin_output, central_species_list=central_species_list,
                        superimpose=superimpose, save_states=save_states, read_states=read_states,
                        diffusion_limited=dflag, check_duplicates=check_duplicates)


if __name__ == '__main__':
    main()
