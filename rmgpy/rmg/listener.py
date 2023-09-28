#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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

import csv
import os

from rmgpy.chemkin import get_species_identifier
from rmgpy.tools.plot import SimulationPlot


class SimulationProfileWriter(object):
    """
    SimulationProfileWriter listens to a ReactionSystem subject
    and writes the species mole numbers as a function of the reaction time
    to a csv file.


    A new instance of the class can be appended to a subject as follows:
    
    reaction_system = ...
    listener = SimulationProfileWriter()
    reaction_system.attach(listener)

    Whenever the subject calls the .notify() method, the
    .update() method of the listener will be called.

    To stop listening to the subject, the class can be detached
    from its subject:

    reaction_system.detach(listener)

    """

    def __init__(self, output_directory, reaction_sys_index, core_species):
        super(SimulationProfileWriter, self).__init__()

        self.output_directory = output_directory
        self.reaction_sys_index = reaction_sys_index
        self.core_species = core_species

    def update(self, reaction_system):
        """
        Opens a file with filename referring to:
            - reaction system
            - number of core species

        Writes to a csv file:
            - header row with species names
            - each row with number of moles of the core species in the given reaction system.
        """

        filename = os.path.join(
            self.output_directory,
            'solver',
            'simulation_{0}_{1:d}.csv'.format(
                self.reaction_sys_index + 1, len(self.core_species)
            )
        )

        header = ['Time (s)', 'Volume (m^3)']
        for spc in self.core_species:
            header.append(get_species_identifier(spc))

        with open(filename, 'w') as csvfile:
            worksheet = csv.writer(csvfile)

            # add header row:
            worksheet.writerow(header)

            # add mole fractions:
            worksheet.writerows(reaction_system.snapshots)


class SimulationProfilePlotter(object):
    """
    SimulationProfilePlotter listens to a ReactionSystem subject
    and plots the top 10 species mole fraction profiles.

    A new instance of the class can be appended to a subject as follows:
    
    reaction_system = ...
    listener = SimulationProfilePlotter()
    reaction_system.attach(listener)

    Whenever the subject calls the .notify() method, the
    .update() method of the listener will be called.

    To stop listening to the subject, the class can be detached
    from its subject:

    reaction_system.detach(listener)
    """

    def __init__(self, output_directory, reaction_sys_index, core_species):
        super(SimulationProfilePlotter, self).__init__()

        self.output_directory = output_directory
        self.reaction_sys_index = reaction_sys_index
        self.core_species = core_species

    def update(self, reaction_system):
        """
        Saves a png with filename referring to:
            - reaction system
            - number of core species
        """

        csv_file = os.path.join(
            self.output_directory,
            'solver',
            'simulation_{0}_{1:d}.csv'.format(
                self.reaction_sys_index + 1, len(self.core_species)
            )
        )

        png_file = os.path.join(
            self.output_directory,
            'solver',
            'simulation_{0}_{1:d}.png'.format(
                self.reaction_sys_index + 1, len(self.core_species)
            )
        )

        SimulationPlot(csv_file=csv_file, num_species=10, ylabel='Moles').plot(png_file)
