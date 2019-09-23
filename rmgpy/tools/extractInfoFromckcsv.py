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
Assume ckcsv contains only one Soln
"""

import csv

import numpy as np

from rmgpy.chemkin import get_species_identifier


def getROPFromCKCSV(ckcsvFile):
    """
    from ckcsv file get three dicts: first_col_dict (e.g., Time/Distance_units:numpy.array)
    spc_total_dict (e.g., species_string: (header,numpy.array)), where header is sth like
    'spc ROP GasRxn Total_units', and 
    spc_indiv_dict (e.g., species_string: list of (header,numpy.array)), where header is 
    sth like 'spc ROP GasRxn#123_units'
    """
    first_col_dict = {}  # store Time/Distance series
    spc_total_dict = {}
    spc_indiv_dict = {}
    read_flag = False
    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            if row[1].strip().startswith('rate-of-production'):
                read_flag = True
            if not read_flag:
                continue
            label = row[0].strip()
            tokens = label.split('_')
            if tokens[0] == 'Time':
                if 'Soln' in tokens[-1]:
                    raise Exception("This function only supports ckcsv with one Soln!")
                units = row[1].strip()[1:-1].lower()
                header = tokens[0] + '_(' + units + ')'

                content_col = np.array([float(r) for r in row[2:]], np.float)
                content_col *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}[units]
                first_col_dict[header] = content_col
                continue

            if tokens[0] == 'Distance':
                if 'Soln' in tokens[-1]:
                    raise Exception("This function only supports ckcsv with one Soln!")

                units = row[1].strip()[1:-1].lower()
                header = tokens[0] + '_(' + units + ')'

                content_col = np.array([float(r) for r in row[2:]], np.float)
                content_col *= {'cm': 1.0, 'mm': 0.1, 'm': 100.}[units]
                first_col_dict[header] = content_col
                continue

            if len(tokens) > 1 and 'ROP' in tokens:
                if 'Soln' in tokens[-1]:
                    raise Exception("This function only supports ckcsv with one Soln!")
                species_string = label.split('_ROP_')[0]
                units = row[1].strip()[1:-1].lower()
                header = ''
                content_col = np.array([float(r) for r in row[2:]], np.float)
                if tokens[-1] == 'Total':
                    header += species_string + ' ROP ' + tokens[2] \
                              + ' ' + tokens[-1] + '_(' + units + ')'
                    if species_string not in spc_total_dict:
                        spc_total_dict[species_string] = (header, content_col)
                    else:
                        raise Exception("ckcsv file has two {} which is not in proper format!".format(species_string))
                else:  # where tokens[-1] is something like GasRxn#123
                    header += species_string + ' ROP ' \
                              + tokens[-1] + '_(' + units + ')'
                    if species_string not in spc_indiv_dict:
                        spc_indiv_dict[species_string] = [(header, content_col)]
                    else:
                        spc_indiv_dict[species_string].append((header, content_col))

    return first_col_dict, spc_total_dict, spc_indiv_dict


def getConcentrationDictFromCKCSV(ckcsvFile):
    """
    from ckcsv file get two dicts: first_col_dict (e.g., Time/Distance_units:numpy.array)
    spc_conc_dict (e.g., species_string: numpy.array)
    """
    first_col_dict = {}  # store Time/Distance series
    spc_conc_dict = {}

    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            tokens = label.split('_')
            if tokens[0] == 'Time':
                if 'Soln' in tokens[-1]:
                    raise Exception("This function only supports ckcsv with one Soln!")
                units = row[1].strip()[1:-1].lower()
                header = tokens[0] + '_(' + units + ')'

                content_col = np.array([float(r) for r in row[2:]], np.float)
                content_col *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}[units]
                first_col_dict[header] = content_col
                continue

            if tokens[0] == 'Distance':
                if 'Soln' in tokens[-1]:
                    raise Exception("This function only supports ckcsv with one Soln!")

                units = row[1].strip()[1:-1].lower()
                header = tokens[0] + '_(' + units + ')'

                content_col = np.array([float(r) for r in row[2:]], np.float)
                content_col *= {'cm': 1.0, 'mm': 0.1, 'm': 100.}[units]
                first_col_dict[header] = content_col
                continue

            if tokens[0] == 'Volume':
                if 'Soln' in tokens[-1]:
                    raise Exception("This function only supports ckcsv with one Soln!")

                units = row[1].strip()[1:-1].lower()
                header = tokens[0] + '_(' + units + ')'

                content_col = np.array([float(r) for r in row[2:]], np.float)
                content_col *= {'cm3': 1.0, 'm3': 1.00e+6}[units]
                first_col_dict[header] = content_col
                continue

            if tokens[0] == 'Temperature':
                if 'Soln' in tokens[-1]:
                    raise Exception("This function only supports ckcsv with one Soln!")

                units = row[1].strip()[1:-1].lower()
                header = tokens[0] + '_(' + units + ')'

                content_col = np.array([float(r) for r in row[2:]], np.float)
                first_col_dict[header] = content_col
                continue

            if tokens[0] == 'Pressure':
                if 'Soln' in tokens[-1]:
                    raise Exception("This function only supports ckcsv with one Soln!")

                units = row[1].strip()[1:-1].lower()
                header = tokens[0] + '_(' + units + ')'

                content_col = np.array([float(r) for r in row[2:]], np.float)
                content_col *= {'bar': 1.0, 'atm': 1.01325}[units]
                first_col_dict[header] = content_col
                continue

            # read concentration (mole fraction profile)
            if len(tokens) > 1:
                if tokens[0] == 'Mole' and tokens[1] == 'fraction':
                    if 'Soln' in tokens[-1]:
                        raise Exception("This function only supports ckcsv with one Soln!")
                    species_string = label.split('Mole_fraction_')[1]
                    content_col = np.array([float(r) for r in row[2:]], np.float)
                    header = species_string + ' Mole_fraction'
                    if species_string not in spc_conc_dict:
                        spc_conc_dict[species_string] = content_col
                    else:
                        raise Exception("ckcsv file has two {} which is not in proper format!".format(species_string))

    return first_col_dict, spc_conc_dict


def getFluxGraphEdgesDict(spc_rop_dict, core_reactions):
    graph_edges_dict = {}
    for rxn in core_reactions:
        for pair in rxn.pairs:
            if pair in graph_edges_dict:
                # get flux from spc_rop_dict
                species_string = get_species_identifier(pair[0])
                flux = getROPFlux(spc_rop_dict, species_string, rxn.index)
                if len(flux) > 0:
                    graph_edges_dict[pair][rxn] = flux
                else:
                    # for rxns like PDD + rad4 == PDD + rad1
                    species_string = get_species_identifier(pair[1])
                    flux = getROPFlux(spc_rop_dict, species_string, rxn.index)
                    if len(flux) > 0:
                        graph_edges_dict[pair][rxn] = -flux
            elif (pair[1], pair[0]) in graph_edges_dict:
                # get flux from spc_rop_dict
                species_string = get_species_identifier(pair[1])
                flux = getROPFlux(spc_rop_dict, species_string, rxn.index)
                if len(flux) > 0:
                    graph_edges_dict[(pair[1], pair[0])][rxn] = flux
                else:
                    species_string = get_species_identifier(pair[0])
                    flux = getROPFlux(spc_rop_dict, species_string, rxn.index)
                    if len(flux) > 0:
                        graph_edges_dict[(pair[1], pair[0])][rxn] = -flux
            else:
                # get flux from spc_rop_dict
                graph_edges_dict[pair] = {}
                species_string = get_species_identifier(pair[0])
                flux = getROPFlux(spc_rop_dict, species_string, rxn.index)
                if len(flux) > 0:
                    graph_edges_dict[pair][rxn] = flux
                else:
                    species_string = get_species_identifier(pair[1])
                    flux = getROPFlux(spc_rop_dict, species_string, rxn.index)
                    if len(flux) > 0:
                        graph_edges_dict[pair][rxn] = -flux
    return graph_edges_dict


def getROPFlux(spc_rop_dict, species_string, rxn_index):
    """
    get the flux (numpy:array) for a given species and given rxn
    """
    if species_string in spc_rop_dict:
        flux_tup_list = spc_rop_dict[species_string]
        for flux_tup in flux_tup_list:
            header = flux_tup[0]
            rxn_num = int(header.split("Rxn#")[1].split('_')[0])
            if rxn_num == rxn_index:
                flux = flux_tup[1]
                return flux
    return []
