#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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

from base64 import b64encode
from io import BytesIO

import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display, HTML


# HTML formatting for output
full = 12
half = full / 2


def generate_header_html(n, fam_rxn, lib_rxn, library_name, families):
    """
    Generates initial lines of HTML for results table.
    """
    html = ['<table style="width:100%;table-layout:fixed;"><tr>']

    if n == 1:
        html += ['<th colspan="{0}" style="color:green">One RMG match for this reaction</th>'.format(full)]
    elif n == 0:
        if families == 'all':
            html += [
                '<th colspan="{0}" style="color:red">Sad :( There are no matches. This is a magic reaction or has chemistry that should be made into a new reaction family</th>'.format(
                    full)]
        else:
            html += [
                '<th colspan="{0}" style="color:red">There are no matches within the selected families: {1}</th>'.format(
                    full, families)]
    else:
        html += [
            '<th colspan="{0}" style="color:blue">There are multiple RMG matches for this reaction. You have to manually create this training reaction.</th>'.format(
                full)]

    html += ['</tr><tr>']
    html += ['<th colspan="{0}">Source Library: {1}</th>'.format(full, library_name)]
    html += ['</tr><tr>']
    if fam_rxn is not None:
        html += ['<td colspan="{0}"><img src="data:image/png;base64,{1}"></td>'.format(full,
                                                                                       b64encode(fam_rxn._repr_png_()))]
    else:
        html += ['<td colspan="{0}"><img src="data:image/png;base64,{1}"></td>'.format(full,
                                                                                       b64encode(lib_rxn._repr_png_()))]
    html += ['</tr><tr>']
    html += ['<th colspan="{0}">Reactant SMILES</th>'.format(half)]
    html += ['<td colspan="{0}">{1}</td>'.format(half, ' + '.join(
        [reactant.molecule[0].toSMILES() for reactant in lib_rxn.reactants]))]
    html += ['</tr><tr>']
    html += ['<th colspan="{0}">Product SMILES</th>'.format(half)]
    html += ['<td colspan="{0}">{1}</td>'.format(half, ' + '.join(
        [product.molecule[0].toSMILES() for product in lib_rxn.products]))]
    html += ['</tr>']

    return html


def generate_template_html(rxn, template):
    """
    Generates HTML for displaying reaction template.
    """
    templateSize = len(template)
    # HTML table uses a 12 column setup, so templates with 5 groups will break the table, which should not happen
    assert templateSize < 5

    html = ['<tr>']
    html += ['<th colspan="{0}">Matched Family</th>'.format(half)]
    html += ['<td colspan="{0}">{1}</td>'.format(half, rxn.family)]
    html += ['</tr><tr>']
    html += ['<th colspan="{0}">Matched Template</th>'.format(half)]
    html += ['<td colspan="{0}">{1}</td>'.format(half, [entry.label for entry in template])]
    html += ['</tr><tr>']
    for entry in template:
        html += ['<td colspan="{0}">{1}</td>'.format(full / templateSize, entry.label)]
    html += ['</tr><tr>']
    for entry in template:
        html += ['<td colspan="{0}"><img src="data:image/png;base64,{1}"></td>'.format(full / templateSize, b64encode(entry.item._repr_png_()))]
    html += ['</tr><tr>']
    if templateSize == 3:
        merged_group = template[0].item.mergeGroups(template[1].item)
        merged_group = merged_group.mergeGroups(template[2].item)
        html += ['<td colspan="{0}">{1}</td>'.format(full, 'Merged Template')]
        html += ['</tr><tr>']
        html += ['<td colspan="{0}"><img src="data:image/png;base64,{1}"></td>'.format(full, b64encode(merged_group._repr_png_()))]

    return html


def process_reactions(database, libraries, families, compareKinetics=True, showAll=False, filterAromatic=True):
    """
    Main function to recreate library reactions from families and display the results.
    """
    master_dict = {}
    multiple_dict = {}

    for library_name in libraries:
        library = database.kinetics.libraries[library_name]
        reaction_dict = {}
        for index, entry in library.entries.iteritems():
            lib_rxn = entry.item
            lib_rxn.kinetics = entry.data
            lib_rxn.index = index

            # Let's make RMG try to generate this reaction from the families!
            fam_rxn_list = database.kinetics.generate_reactions_from_families(
                reactants=lib_rxn.reactants,
                products=lib_rxn.products,
                only_families=None if families == 'all' else families,
                resonance=True,
            )

            # Filter by aromatic resonance structures if requested
            if filterAromatic and len(fam_rxn_list) > 1:
                selected_rxns = []
                max_num_aromatic_reactants = 0
                for fam_rxn in fam_rxn_list:
                    num_aromatic_reactants = 0
                    reactants = fam_rxn.reactants if fam_rxn.isForward else fam_rxn.products
                    for r in reactants:
                        num_aromatic_reactants += r.molecule[0].isAromatic()
                    if num_aromatic_reactants > max_num_aromatic_reactants:
                        max_num_aromatic_reactants = num_aromatic_reactants
                        selected_rxns = [fam_rxn]
                    elif num_aromatic_reactants == max_num_aromatic_reactants:
                        selected_rxns.append(fam_rxn)
                    else:
                        continue
                if selected_rxns:
                    fam_rxn_list = selected_rxns

            if len(fam_rxn_list) == 1:
                fam_rxn = fam_rxn_list[0]

                forward = fam_rxn.isForward

                # Find the labeled atoms using family and reactants & products from fam_rxn
                database.kinetics.families[fam_rxn.family].addAtomLabelsForReaction(fam_rxn)

                # Replace lib_rxn spcs with fam_rxn spcs to transfer atom labels
                if forward:
                    lib_rxn.reactants = fam_rxn.reactants
                    lib_rxn.products = fam_rxn.products
                    lib_rxn._degeneracy = fam_rxn.degeneracy
                else:
                    lib_rxn.reactants = fam_rxn.products
                    lib_rxn.products = fam_rxn.reactants

                if fam_rxn.family in reaction_dict:
                    reaction_dict[fam_rxn.family].append(lib_rxn)
                else:
                    reaction_dict[fam_rxn.family] = [lib_rxn]

                template = database.kinetics.families[fam_rxn.family].retrieveTemplate(fam_rxn.template)

                if compareKinetics:
                    # Check what the current kinetics for this template are
                    newKinetics = lib_rxn.kinetics
                    oldKinetics = database.kinetics.families[fam_rxn.family].getKineticsForTemplate(template, degeneracy=fam_rxn.degeneracy)[0]
                    # Evaluate kinetics
                    tlistinv = np.linspace(1000 / 1500, 1000 / 300, num=10)
                    tlist = 1000 * np.reciprocal(tlistinv)
                    newklist = np.log10(np.array([newKinetics.getRateCoefficient(t) for t in tlist]))
                    oldklist = np.log10(np.array([oldKinetics.getRateCoefficient(t) for t in tlist]))
                    # Create plot
                    plt.cla()
                    plt.plot(tlistinv, newklist, label='New')
                    plt.plot(tlistinv, oldklist, label='Current')
                    plt.legend()
                    plt.xlabel('1000/T')
                    plt.ylabel('log(k)')
                    fig = BytesIO()
                    plt.savefig(fig, format='png')
                    fig.seek(0)
                    figdata = b64encode(fig.getvalue())
                    fig.close()

                # Format output using html
                html = generate_header_html(1, fam_rxn, lib_rxn, library_name, families)
                html += generate_template_html(fam_rxn, template)
                if compareKinetics:
                    if not forward:
                        html += ['<tr><th colspan="{0}" style="color:red;text-align:center">'
                                 'Note: Training reaction written in opposite direction from reaction family.'
                                 '</th></tr>'.format(full)]
                    html += ['<tr>']
                    html += ['<td colspan="{0}"><strong>New Kinetics:</strong><br>{1}<br><br>'
                             '<strong>Current Kinetics</strong><br>{2}</td>'.format(half, newKinetics, oldKinetics)]
                    html += ['<td colspan="{0}"><img src="data:image/png;base64,{1}"></td>'.format(half, figdata)]
                    html += ['</tr>']
                html += ['</table>']

                display(HTML(''.join(html)))
            elif len(fam_rxn_list) == 0:
                if showAll:
                    html = generate_header_html(0, None, lib_rxn, library_name, families)
                    html += ['</table>']

                    display(HTML(''.join(html)))
                else:
                    continue
            else:
                # Save results to allow further processing later
                if library_name in multiple_dict:
                    multiple_dict[library_name].append((lib_rxn, fam_rxn_list))
                else:
                    multiple_dict[library_name] = [(lib_rxn, fam_rxn_list)]

                if compareKinetics: oldKinetics = []

                for i, rxn in enumerate(fam_rxn_list):
                    forward = rxn.isForward

                    template = database.kinetics.families[rxn.family].retrieveTemplate(rxn.template)

                    if compareKinetics:
                        oldKinetics.append(database.kinetics.families[rxn.family].getKineticsForTemplate(template, degeneracy=rxn.degeneracy)[0])

                    if i == 0:
                        html = generate_header_html(2, rxn, lib_rxn, library_name, families)

                    html += ['<tr>']
                    html += ['<th colspan="{0}">Match #{1} - For the following resonance form of the reaction:</th>'.format(full, i + 1)]
                    html += ['</tr><tr>']
                    html += ['<td colspan="{0}"><img src="data:image/png;base64,{1}"></td>'.format(full, b64encode(rxn._repr_png_()))]
                    html += ['</tr>']
                    html += generate_template_html(rxn, template)

                if compareKinetics:
                    newKinetics = lib_rxn.kinetics
                    # Evaluate kinetics
                    tlistinv = np.linspace(1000 / 1500, 1000 / 300, num=10)
                    tlist = 1000 * np.reciprocal(tlistinv)
                    newklist = np.log10(np.array([newKinetics.getRateCoefficient(t) for t in tlist]))
                    oldklist = []
                    for kinetics in oldKinetics:
                        oldklist.append(np.log10(np.array([kinetics.getRateCoefficient(t) for t in tlist])))
                    # Create plot
                    plt.cla()
                    plt.plot(tlistinv, newklist, label='New')
                    for i, k in enumerate(oldklist):
                        plt.plot(tlistinv, k, label='Match #{0}'.format(i + 1))
                    plt.legend()
                    plt.xlabel('1000/T')
                    plt.ylabel('log(k)')
                    fig = BytesIO()
                    plt.savefig(fig, format='png')
                    fig.seek(0)
                    figdata = b64encode(fig.getvalue())
                    fig.close()

                    if not forward:
                        html += ['<tr><th colspan="{0}" style="color:red;text-align:center">'
                                 'Note: Training reaction written in opposite direction from reaction family.'
                                 '</tr></tr>'.format(full)]
                    html += ['<tr><td colspan="{0}">'.format(half)]
                    html += ['<strong>New Kinetics:</strong><br>{0}'.format(newKinetics)]
                    for i, kinetics in enumerate(oldKinetics):
                        html += ['<br><br><strong>Match #{0} Kinetics:</strong><br>{1}'.format(i + 1, kinetics)]
                    html += ['</td><td colspan="{0}"><img src="data:image/png;base64,{1}"></td>'.format(half, figdata)]
                    html += ['</tr>']

                html += ['</table>']

                display(HTML(''.join(html)))

        # Save results for this library
        if reaction_dict:
            master_dict[library_name] = reaction_dict

    return master_dict, multiple_dict


def review_reactions(master_dict, prompt=False):
    """
    Function to display reactions that will be added to training depositories.

    If prompt is True, will ask the user whether each reaction shoule be added.
    """

    print('================================================================================')
    print('The following reactions will be added to the indicated families.')

    for library_name, reaction_dict in master_dict.iteritems():
        print('================================================================================')
        print('Source Library: {0}'.format(library_name))
        for family_name, reaction_list in reaction_dict.iteritems():
            print('--------------------------------------------------------------------------------')
            print('Destination Family: {0}'.format(family_name))

            index = 0
            while index < len(reaction_list):
                reaction = reaction_list[index]
                print('\nOriginal Library Reaction: {0}'.format(reaction.kinetics.comment))
                display(reaction)
                print(reaction.kinetics)

                if prompt:
                    success = False
                    while not success:
                        choice = raw_input('Would you like to add this reaction? (y/n) ')
                        if choice in ['y', 'n']:
                            success = True
                        else:
                            print('Invalid choice.')

                    if choice == 'y':
                        index += 1
                    elif choice == 'n':
                        del reaction_list[index]
                else:
                    index += 1

    print('================================================================================')
    print('All reactions reviewed.')
    print('================================================================================')


def manual_selection(master_dict, multiple_dict, database):
    """
    For reactions with multiple matches, prompts the user to choose one match
    to be added to the training depositories.
    """

    print('================================================================================')
    print('The following reactions had multiple matches. You may choose one match to add.')

    for library_name, reaction_list in multiple_dict.iteritems():
        print('================================================================================')
        print('Source Library: {0}'.format(library_name))

        index = 0
        while index < len(reaction_list):
            lib_rxn, fam_rxn_list = reaction_list[index]
            print('--------------------------------------------------------------------------------')
            print('Original Library Reaction: {0}\n'.format(lib_rxn.kinetics.comment))

            for i, fam_rxn in enumerate(fam_rxn_list):
                print('Match #{0}'.format(i))
                print('Reaction Family: {0}'.format(fam_rxn.family))
                print('Reaction Template: {0}'.format(fam_rxn.template))
                display(fam_rxn)

            success = False
            while not success:
                choice = raw_input('Select a match to add (or use "s" to skip): ')
                try:
                    choice = int(choice)
                except ValueError:
                    pass
                if choice in ['s', 'q'] or (choice >= 0 and choice <= i):
                    success = True
                else:
                    print('Invalid choice.')

            if choice == 's':
                print('Skipping this reaction.')
                index += 1
                continue

            print('Adding match #{0} to list of new training reactions.'.format(choice))

            fam_rxn = fam_rxn_list[choice]

            forward = fam_rxn.isForward

            # Find the labeled atoms using family and reactants & products from fam_rxn
            database.kinetics.families[fam_rxn.family].addAtomLabelsForReaction(fam_rxn)

            # Replace lib_rxn spcs with fam_rxn spcs to transfer atom labels
            if forward:
                lib_rxn.reactants = fam_rxn.reactants
                lib_rxn.products = fam_rxn.products
                lib_rxn._degeneracy = fam_rxn.degeneracy
            else:
                lib_rxn.reactants = fam_rxn.products
                lib_rxn.products = fam_rxn.reactants

            try:
                reaction_dict = master_dict[library_name]
            except KeyError:
                reaction_dict = {}
                master_dict[library_name] = reaction_dict

            if fam_rxn.family in reaction_dict:
                reaction_dict[fam_rxn.family].append(lib_rxn)
            else:
                reaction_dict[fam_rxn.family] = [lib_rxn]

            # Remove this item from the list to prevent reprocessing
            del reaction_list[index]

    print('================================================================================')
    print('Manual selection of reactions completed.')
    print('================================================================================')

