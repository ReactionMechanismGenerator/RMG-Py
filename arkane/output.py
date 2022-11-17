#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains helper functionality for writing Arkane output files.
"""

import ast
import logging
import os
import shutil

from rmgpy.data.base import Entry
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.species import Species

################################################################################


class PrettifyVisitor(ast.NodeVisitor):
    """
    A class for traversing an abstract syntax tree to assemble a prettier
    version of the code used to create the tree. Used by the :func:`prettify`
    function.
    """

    def __init__(self, level=0, indent=4):
        self.string = ''
        self.level = level
        self.indent = indent

    def visit_Call(self, node):
        """
        Return a pretty representation of the class or function call represented by `node`.
        """
        keywords = []
        for keyword in node.keywords:
            keywords.append('{0}={1}'.format(keyword.arg, self.visit(keyword.value)))
        result = '{0}({1})'.format(node.func.id, ', '.join(keywords))

        if len(result) > 80:
            result = node.func.id + '(\n'

            self.level += 1
            for keyword in node.keywords:
                result += '{2}{0} = {1},\n'.format(keyword.arg, self.visit(keyword.value),
                                                   ' ' * (self.level * self.indent))
            self.level -= 1

            result += ' ' * (self.level * self.indent) + ')'

        self.string = result

        return result

    def visit_List(self, node):
        """
        Return a pretty representation of the list represented by `node`.
        """
        if any([not isinstance(e, (ast.Str, ast.Num, ast.UnaryOp)) for e in node.elts]):
            # Split elements onto multiple lines
            result = '[\n'
            self.level += 1
            for element in node.elts:
                result += '{1}{0},\n'.format(self.visit(element), ' ' * (self.level * self.indent))
            self.level -= 1
            result += '{0}]'.format(' ' * (self.level * self.indent))
            return result
        else:
            # Keep elements on one line
            result = '[{0}]'.format(', '.join([self.visit(e) for e in node.elts]))
            self.string = result
            return result

    def visit_Tuple(self, node):
        """
        Return a pretty representation of the tuple represented by `node`.
        """
        # If the tuple represents a quantity, keep it on one line
        is_quantity = True
        if len(node.elts) == 0 or not isinstance(node.elts[0], (ast.Num, ast.List)) or (
                isinstance(node.elts[0], ast.List) and
                any([not isinstance(e, (ast.Num, ast.UnaryOp)) for e in node.elts[0].elts])):
            is_quantity = False
        elif len(node.elts) < 2 or not isinstance(node.elts[1], ast.Str):
            is_quantity = False

        if not is_quantity:
            # Split elements onto multiple lines
            result = '(\n'
            self.level += 1
            for element in node.elts:
                result += '{1}{0},\n'.format(self.visit(element), ' ' * (self.level * self.indent))
            self.level -= 1
            result += '{0})'.format(' ' * (self.level * self.indent))
            return result
        else:
            # Keep elements on one line
            result = '({0})'.format(', '.join([self.visit(e) for e in node.elts]))
            self.string = result
            return result

    def visit_Dict(self, node):
        """
        Return a pretty representation of the dict represented by `node`.
        """
        if (any([not isinstance(e, (ast.Str, ast.Num)) for e in node.keys])
                or any([not isinstance(e, (ast.Str, ast.Num)) for e in node.values])):
            # Split elements onto multiple lines
            result = '{\n'
            self.level += 1
            for key, value in zip(node.keys, node.values):
                result += '{2}{0}: {1},\n'.format(self.visit(key), self.visit(value), ' ' * (self.level * self.indent))
            self.level -= 1
            result += '{0}}}'.format(' ' * (self.level * self.indent))
            self.string = result
            return result
        else:
            # Keep elements on one line
            result = '{{{0}}}'.format(', '.join(['{0}: {1}'.format(self.visit(key), self.visit(value))
                                                 for key, value in zip(node.keys, node.values)]))
            self.string = result
            return result

    def visit_Str(self, node):
        """
        Return a pretty representation of the string represented by `node`.
        """
        result = repr(node.s)
        self.string = result
        return result

    def visit_Num(self, node):
        """
        Return a pretty representation of the number represented by `node`.
        """
        result = '{0:g}'.format(node.n)
        # result = repr(node.n)
        self.string = result
        return result

    def visit_UnaryOp(self, node):
        """
        Return a pretty representation of the number represented by `node`.
        """
        operators = {
            ast.UAdd: '+',
            ast.USub: '-',
        }
        result = '{0}{1}'.format(operators[node.op.__class__], self.visit(node.operand))
        self.string = result
        return result


def prettify(string, indent=4):
    """
    Return a "pretty" version of the given `string`, representing a snippet of
    Python code such as a representation of an object or function. This 
    involves splitting of tuples, lists, and dicts (including parameter lists)
    onto multiple lines, indenting as appropriate for readability.
    """
    # Parse the node into an abstract syntax tree
    node = ast.parse(string)
    # Traverse the tree, assembling the pretty version of the string
    visitor = PrettifyVisitor(indent=indent)
    visitor.visit(node)
    # Return the pretty version of the string
    return visitor.string


def get_str_xyz(spc):
    """
    Get a string representation of the 3D coordinates from the conformer.

    Args:
        spc (Species): A Species instance.

    Returns:
        str: A string representation of the coordinates
    """
    if spc.conformer.coordinates is not None:
        from arkane.common import symbol_by_number
        xyz_list = list()
        for number, coord in zip(spc.conformer.number.value_si, spc.conformer.coordinates.value_si):
            coord_angstroms = coord * 10 ** 10
            row = f'{symbol_by_number[number]:4}'
            row += '{0:14.8f}{1:14.8f}{2:14.8f}'.format(*coord_angstroms)
            xyz_list.append(row)
        return '\n'.join(xyz_list)
    else:
        return None


def save_thermo_lib(species_list, path, name, lib_long_desc):
    """
    Save an RMG thermo library.

    Args:
        species_list (list): Entries are Species object instances for which thermo will be saved.
        path (str): The base folder in which the thermo library will be saved.
        name (str): The library name.
        lib_long_desc (str): A multiline string with relevant description.
    """
    if species_list:
        lib_path = os.path.join(path, f'{name}.py')
        thermo_library = ThermoLibrary(name=name, long_desc=lib_long_desc)
        for i, spc in enumerate(species_list):
            if spc.thermo is not None:
                long_thermo_description = f'\nSpin multiplicity: {spc.conformer.spin_multiplicity}' \
                                          f'\nExternal symmetry: {spc.molecule[0].symmetry_number}' \
                                          f'\nOptical isomers: {spc.conformer.optical_isomers}\n'
                xyz = get_str_xyz(spc)
                if xyz is not None:
                    long_thermo_description += f'\nGeometry:\n{xyz}'
                thermo_library.load_entry(index=i,
                                          label=spc.label,
                                          molecule=spc.molecule[0].to_adjacency_list(),
                                          thermo=spc.thermo,
                                          shortDesc=spc.thermo.comment,
                                          longDesc=long_thermo_description)
            else:
                logging.warning(f'Species {spc.label} did not contain any thermo data and was omitted from the thermo '
                                f'library {name}.')
        thermo_library.save(lib_path)


def save_kinetics_lib(rxn_list, path, name, lib_long_desc):
    """
    Save an RMG kinetics library.

    Args:
        rxn_list (list): Entries are Reaction object instances for which kinetics will be saved.
        path (str): The base folder in which the kinetic library will be saved.
        name (str): The library name.
        lib_long_desc (str): A multiline string with relevant description.
    """
    entries = dict()
    if rxn_list:
        for i, rxn in enumerate(rxn_list):
            if rxn.kinetics is not None:
                entry = Entry(
                    index=i,
                    item=rxn,
                    data=rxn.kinetics,
                    label=' <=> '.join([' + '.join([reactant.label for reactant in rxn.reactants]),
                                        ' + '.join([product.label for product in rxn.products])]),
                )
                entries[i+1] = entry
            else:
                logging.warning(f'Reaction {rxn.label} did not contain any kinetic data and was omitted from the '
                                f'kinetics library.')
        kinetics_library = KineticsLibrary(name=name, long_desc=lib_long_desc, auto_generated=True)
        kinetics_library.entries = entries
        if os.path.exists(path):
            shutil.rmtree(path)
        try:
            os.makedirs(path)
        except OSError:
            pass
        kinetics_library.save(os.path.join(path, 'reactions.py'))
        kinetics_library.save_dictionary(os.path.join(path, 'dictionary.txt'))
