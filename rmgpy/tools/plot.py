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

from __future__ import division

import os

import matplotlib as mpl

# Force matplotlib to not use any Xwindows backend.
# This must be called before pylab, matplotlib.pyplot, or matplotlib.backends is imported
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from rmgpy.tools.data import GenericData


def plot_sensitivity(outputDirectory, reactionSystemIndex, sensitiveSpeciesList, number=10, fileformat='.png'):
    """
    A function for plotting the top reaction thermo sensitivities (the number is
    inputted as the variable `number`) in bar plot format.
    To be called after running a simulation on a particular reactionSystem.
    """

    for species in sensitiveSpeciesList:
        csv_file = os.path.join(
            outputDirectory,
            'solver',
            'sensitivity_{0}_SPC_{1}.csv'.format(
                reactionSystemIndex + 1, species.index
            )
        )

        reaction_plot_file = os.path.join(
            outputDirectory,
            'solver',
            'sensitivity_{0}_SPC_{1}_reactions'.format(
                reactionSystemIndex + 1, species.index
            ) + fileformat
        )

        thermo_plot_file = os.path.join(
            outputDirectory,
            'solver',
            'sensitivity_{0}_SPC_{1}_thermo'.format(
                reactionSystemIndex + 1, species.index
            ) + fileformat
        )

        ReactionSensitivityPlot(csvFile=csv_file, numReactions=number).barplot(reaction_plot_file)
        ThermoSensitivityPlot(csvFile=csv_file, numSpecies=number).barplot(thermo_plot_file)


def parseCSVData(csvFile):
    """
    This function parses a typical csv file outputted from a simulation or
    sensitivity analysis in the form of
    
    Time (s)  Header1  Header2  Header....
    t0        val1_0   val2_0   val...
    t1..
    ..
    
    It returns the data in the form 
    
    Time, DataList
    
    Where Time is returned as a GenericData object, and DataList is list of GenericData objects
    """
    import csv
    import re

    # Pattern for matching indices or units
    index_pattern = re.compile(r'^\S+\(\d+\)$')
    units_pattern = re.compile(r'\s\(.+\)$')
    rxn_sens_pattern = re.compile('^dln\[\S+\]\/dln\[k\d+\]:\s\S+$')
    thermo_sens_pattern = re.compile('^dln\[\S+\]\/dG\[\S+\]$')

    time_data = []
    data = {}
    f = csv.reader(open(csvFile, 'r'))

    columns = list(zip(*f))
    time = GenericData(label=columns[0][0],
                       data=np.array(columns[0][1:], dtype=np.float64),
                       )

    # Parse the units from the Time header
    if units_pattern.search(time.label):
        label, sep, units = time.label[:-1].rpartition('(')
        time.label = label
        time.units = units

    data_list = []
    for col in columns[1:]:
        header = col[0]
        values = np.array(col[1:], dtype=np.float64)
        data = GenericData(label=header, data=values)

        # Parse the index or the label from the header
        if index_pattern.search(data.label):
            species, sep, index = data.label[:-1].rpartition('(')
            # Save the species attribute if an index was found
            data.species = species
            data.index = int(index)
        elif units_pattern.search(data.label):
            label, sep, units = data.label[:-1].rpartition('(')
            data.label = label
            data.units = units
        elif rxn_sens_pattern.search(data.label):
            rxn = data.label.split()[1]
            index = data.label.split()[0][:-2].rpartition('dln[k')[2]
            data.reaction = rxn
            data.index = int(index)
        elif thermo_sens_pattern.search(data.label):
            species = data.label[:-1].rpartition('dG[')[2]
            data.species = species
            if index_pattern.search(species):
                data.index = int(species[:-1].rpartition('(')[2])

        data_list.append(data)

    return time, data_list


def findNearest(array, value):
    """
    Returns the index of the closest value in a sorted array
    """
    idx = (np.abs(array - value)).argmin()
    return idx


def linearlyInterpolatePoint(xArray, yArray, xValue):
    """
    Returns the interpolated yValue for given xValue using data from the two sorted arrays:
    """
    # Find the next largest point in xArray that is still smaller than xValue:
    lower_index = None
    for index, x in enumerate(xArray):
        if x > xValue:
            break
        lower_index = index

    # If xValue is outside the domain of xArray, we use either the min or max points for dydx
    if lower_index is None:
        lower_index = 0
    elif lower_index == len(xArray) - 1:
        lower_index = lower_index - 1
    higher_index = lower_index + 1

    dydx = (yArray[higher_index] - yArray[lower_index]) / (xArray[higher_index] - xArray[lower_index])

    if xValue < xArray[lower_index]:
        yValue = yArray[lower_index] - dydx * (xValue - xArray[lower_index])
    else:
        yValue = yArray[lower_index] + dydx * (xValue - xArray[lower_index])
    return yValue


class GenericPlot(object):
    """
    A generic plotting class that can be extended to plot other things.
    """

    def __init__(self, xVar=None, yVar=None, title='', xlabel='', ylabel=''):
        self.xVar = xVar
        # Convert yVar to a list if it wasn't one already
        if isinstance(yVar, GenericData):
            self.yVar = [yVar]
        else:
            self.yVar = yVar
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel

    def plot(self, filename=''):
        """
        Execute the actual plotting
        """
        mpl.rc('font', family='sans-serif')
        fig = plt.figure()

        ax = fig.add_subplot(111)

        x_var = self.xVar
        y_var = self.yVar

        if len(y_var) == 1:
            y = y_var[0]
            ax.plot(x_var.data, y.data)

            # Create a ylabel for the label of the y variable
            if not self.ylabel and y.label:
                ylabel = y.label
                if y.units: ylabel += ' ({0})'.format(y.units)
                plt.ylabel(ylabel)
        else:
            for y in y_var:
                ax.plot(x_var.data, y.data, '-', label=y.label)

        if self.xlabel:
            plt.xlabel(self.xlabel)
        elif x_var.label:
            xlabel = x_var.label
            if x_var.units: xlabel += ' ({0})'.format(x_var.units)
            plt.xlabel(xlabel)

        if self.ylabel:
            plt.ylabel(self.ylabel)

        if self.title:
            plt.title(self.title)

        ax.grid(True)
        handles, labels = ax.get_legend_handles_labels()
        if labels:
            # Create a legend outside the plot and adjust width based off of longest legend label
            max_string_length = max([len(label) for label in labels])
            width = 1.05 + .011 * max_string_length
            legend = ax.legend(handles, labels, loc='upper center', numpoints=1,
                               bbox_to_anchor=(width, 1))  # bbox_to_anchor=(1.01,.9)
            fig.savefig(filename, bbox_extra_artists=(legend,), bbox_inches='tight')
        else:
            fig.savefig(filename, bbox_inches='tight')

    def barplot(self, filename='', idx=None):
        """
        Plot a generic barplot using just the yVars.
        idx is the index of the each y-variable to be plotted. if not given, the last value will be used
        """
        mpl.rc('font', family='sans-serif')

        fig = plt.figure()
        ax = fig.add_subplot(111)

        position = np.arange(len(self.yVar), 0, -1)
        # Reverse in order to go front top to bottom
        if not idx:
            idx = -1
        ax.barh(position, np.array([y.data[idx] for y in self.yVar]), align='center', alpha=0.5)
        plt.yticks(position, [y.label for y in self.yVar])

        # If any labels or titles are explicitly specified, write them
        if self.xlabel:
            plt.xlabel(self.xlabel)

        if self.ylabel:
            plt.ylabel(self.ylabel)

        if self.title:
            plt.title(self.title)

        plt.axis('tight')
        fig.savefig(filename, bbox_inches='tight')

    def comparePlot(self, otherGenericPlot, filename='', title='', xlabel='', ylabel=''):
        """
        Plot a comparison data plot of this data vs a second GenericPlot class
        """

        mpl.rc('font', family='sans-serif')
        # mpl.rc('text', usetex=True)
        fig = plt.figure()

        ax = fig.add_subplot(111)

        styles = ['-', ':']
        # Plot the sets of data
        for i, plot in enumerate([self, otherGenericPlot]):
            # Reset the color cycle per plot to get matching colors in each set
            plt.gca().set_prop_cycle(None)

            x_var = plot.xVar
            y_var = plot.yVar
            # Convert y_var to a list if it wasn't one already
            if isinstance(y_var, GenericData):
                y_var = [y_var]

            if len(y_var) == 1:
                y = y_var[0]
                ax.plot(x_var.data, y.data, styles[i])
                # Save a ylabel based on the y variable's label if length of y_var contains only 1 variable
                if not self.ylabel and y.label:
                    self.ylabel = y.label
                    if y.units:
                        self.ylabel += ' ({0})'.format(y.units)
            else:
                for y in y_var:
                    ax.plot(x_var.data, y.data, styles[i], label=y.label)

            # Plot the second set of data

        # Prioritize using the function's x and y labels, otherwise the labels from this data object
        if xlabel:
            plt.xlabel(xlabel)
        elif self.xlabel:
            plt.xlabel(self.xlabel)
        elif self.xVar.label:
            xlabel = self.xVar.label
            if self.xVar.units: xlabel += ' ({0})'.format(self.xVar.units)
            plt.xlabel(xlabel)

        if ylabel:
            plt.ylabel(ylabel)
        elif self.ylabel:
            plt.ylabel(self.ylabel)

        # Use user inputted title
        if title:
            plt.title(title)

        ax.grid(True)
        handles, labels = ax.get_legend_handles_labels()
        if labels:
            # Create a legend outside the plot and adjust width based off of longest legend label
            max_string_length = max([len(label) for label in labels])
            width = 1.2 + .011 * max_string_length * 2
            legend = ax.legend(handles, labels, loc='upper center', numpoints=1,
                               bbox_to_anchor=(width, 1), ncol=2)  # bbox_to_anchor=(1.01,.9)
            fig.savefig(filename, bbox_extra_artists=(legend,), bbox_inches='tight')
        else:
            fig.savefig(filename, bbox_inches='tight')


class SimulationPlot(GenericPlot):
    """
    A class for plotting simulations containing mole fraction vs time data. 
    Can plot the top species in generic simulation csv generated by RMG-Py
    i.e. simulation_1_19.csv, found in the solver folder of an RMG job
    
    Use numSpecies as a flag to dictate how many species to plot.  
    This function will plot the top species, based on maximum mole fraction at
    any point in the simulation.
    
    Alternatively, the `species` flag can be used as a dictionary for 
    plotting specific species within the csvFile
    This should be formulated as 
    {'desired_name_for_species': 'corresponding_chemkin_name_of_species'}
    """

    def __init__(self, xVar=None, yVar=None, title='', xlabel='', ylabel='', csvFile='', numSpecies=None, species=None):
        GenericPlot.__init__(self, xVar=xVar, yVar=yVar, title=title, xlabel=xlabel, ylabel=ylabel)
        self.csvFile = csvFile
        self.numSpecies = numSpecies
        self.species = species if species else {}

    def load(self):
        if self.xVar == None and self.yVar == None:
            time, data_list = parseCSVData(self.csvFile)
        else:
            time = self.xVar
            data_list = self.yVar

        species_data = []
        if self.species:
            # A specific set of species was specified to be plotted
            for speciesLabel, chemkinLabel in self.species.items():
                for data in data_list:
                    if chemkinLabel == data.label:
                        # replace the data label with the desired species label
                        data.label = speciesLabel
                        species_data.append(data)
                        break

        else:
            for data in data_list:
                # Only plot if RMG detects that the data corresponds with a species
                # This will not include bath gases
                if data.species:
                    species_data.append(data)

        self.xVar = time
        self.yVar = species_data

    def plot(self, filename=''):
        filename = filename if filename else 'simulation.png'
        self.load()
        self.yVar.sort(key=lambda x: max(x.data), reverse=True)

        if self.numSpecies:
            self.yVar = self.yVar[:self.numSpecies]

        GenericPlot.plot(self, filename=filename)

    def comparePlot(self, otherSimulationPlot, filename='', title='', xlabel='', ylabel=''):

        filename = filename if filename else 'simulation_compare.png'
        self.load()
        otherSimulationPlot.load()

        # Restrict the number of species
        if self.numSpecies:
            self.yVar = self.yVar[:self.numSpecies]
            otherSimulationPlot.yVar = otherSimulationPlot.yVar[:self.numSpecies]
        GenericPlot.comparePlot(self, otherSimulationPlot, filename, title, xlabel, ylabel)


class ReactionSensitivityPlot(GenericPlot):
    """
    A class for plotting the top reaction sensitivites in a generic sensitivity csv file generated by RMG-Py.
    
    `numReactions` is a flag indicating the number of reaction sensitivities to plot.
    This function will plot the top sensitivities based on this number, based on the
    magnitude of the sensitivity at the final time step.
    
    `reactions` can be used to plot the sensitivities of a specific set of reactions.
    This should be formulated as 
    {'desired_rxn_string': 'corresponding chemkin rxn string'}
    
    barplot() will instead plot a horizontal bar plot of the sensitivities at a given
    time step.  If time step is not given, the end step will automatically be chosen
    """

    def __init__(self, xVar=None, yVar=None, title='', xlabel='', ylabel='', csvFile='', numReactions=None,
                 reactions=None):
        GenericPlot.__init__(self, xVar=xVar, yVar=yVar, title=title, xlabel=xlabel, ylabel=ylabel)
        self.csvFile = csvFile
        self.numReactions = numReactions
        self.reactions = reactions if reactions else {}

    def load(self):
        if self.xVar == None and self.yVar == None:
            time, data_list = parseCSVData(self.csvFile)
        else:
            time = self.xVar
            data_list = self.yVar
        reaction_data = []
        if self.reactions:
            # A specific set of reaction sensitivities was specified to be plotted
            for reaction_label, chemkin_label in self.reactions.items():
                for data in data_list:
                    if chemkin_label == data.reaction:
                        # replace the data label with the desired species label
                        data.label = reaction_label
                        reaction_data.append(data)
                        break
        else:
            for data in data_list:
                if data.reaction:
                    reaction_data.append(data)
        self.xVar = time
        self.yVar = reaction_data

    def plot(self, filename=''):
        filename = filename if filename else "reaction_sensitivity.png"
        self.load()
        # Sort reactions according to absolute max value of final time point
        self.yVar.sort(key=lambda x: abs(x.data[-1]), reverse=True)
        if self.numReactions:
            self.yVar = self.yVar[:self.numReactions]
        self.ylabel = 'dln(c)/dln(k_i)'
        GenericPlot.plot(self, filename=filename)

    def barplot(self, filename='', t=None):
        """
        Time must be indicated in seconds
        The closest timepoint will then be used, otherwise if no time point is given,
        the end time step will be used
        """
        filename = filename if filename else "reaction_sensitivity.png"
        self.load()
        # Sort reactions according to absolute max value at the specified time point
        # if the time point is not given, use the final time point
        if t:
            idx = findNearest(self.xVar.data, t)
        else:
            idx = -1

        self.yVar.sort(key=lambda x: abs(x.data[idx]), reverse=True)
        if self.numReactions:
            self.yVar = self.yVar[:self.numReactions]

        if not self.xlabel:
            self.xlabel = 'dln(c)/dln(k_i)'
        GenericPlot.barplot(self, filename=filename, idx=idx)

    def uncertaintyPlot(self, totalVariance, t=None, filename=''):
        """
        Plot the top uncertainty contributions resulting from uncertainties in the
        kinetic parameters.  The totalVariance must be specified.  Optionally,
        the reaction time `t` in seconds can be specified for plotting the uncertainties.
        The number of reaction uncertainties to plot is determined by self.numReactions
        """

        filename = filename if filename else "kinetics_uncertainty.png"
        self.load()
        if t:
            idx = findNearest(self.xVar.data, t)
        else:
            idx = -1

        reaction_uncertainty = []

        total_uncertainty = totalVariance

        for reactionSens in self.yVar:
            if isinstance(reactionSens, np.ndarray):
                # The parameter uncertainties are an array which should have the same length as the sensitivity data
                uncertainty_data = reactionSens * reactionSens.uncertainty
                uncertainty_contribution = uncertainty_data[idx] ** 2
            else:
                uncertainty_contribution = (reactionSens.data[idx] * reactionSens.uncertainty) ** 2

            reaction_uncertainty.append([reactionSens.label, reactionSens.reaction, uncertainty_contribution])

        # Normalize and create new list of GenericData
        new_y_var = []
        for label, reaction, uncertainty in reaction_uncertainty:
            data = GenericData(label=label, reaction=reaction, data=[uncertainty / total_uncertainty * 100])
            new_y_var.append(data)

        new_y_var.sort(key=lambda x: abs(x.data[0]), reverse=True)
        new_y_var = new_y_var[:self.numReactions]

        GenericPlot(xVar=None, yVar=new_y_var, xlabel="Uncertainty Contribution (%)").barplot(filename=filename)

        return reaction_uncertainty


class ThermoSensitivityPlot(GenericPlot):
    """
    A class for plotting the top sensitivities to a thermo DeltaG value of species within the model.
    The value used is the sensitivity at the final time point.
    `numSpecies` indicates the number of species to plot. 
    
    `species` is a dictionary corresponding to specific species thermo sensitivities to be plotted
    
    barplot() will instead plot a horizontal bar plot of the sensitivities at a given
    time step.  If time step is not given, the end step will automatically be chosen
    """

    def __init__(self, xVar=None, yVar=None, title='', xlabel='', ylabel='', csvFile='', numSpecies=None, species=None):
        GenericPlot.__init__(self, xVar=xVar, yVar=yVar, title=title, xlabel=xlabel, ylabel=ylabel)
        self.csvFile = csvFile
        self.numSpecies = numSpecies
        self.species = species if species else {}

    def load(self):
        if self.xVar == None and self.yVar == None:
            time, data_list = parseCSVData(self.csvFile)
        else:
            time = self.xVar
            data_list = self.yVar
        thermo_data = []
        if self.species:
            # A specific set of species sensitivities was specified to be plotted
            for speciesLabel, chemkinLabel in self.species.items():
                for data in data_list:
                    if chemkinLabel == data.species:
                        # replace the data label with the desired species label
                        data.label = speciesLabel
                        thermo_data.append(data)
                        break
        else:
            for data in data_list:
                if data.species:
                    thermo_data.append(data)

        self.xVar = time
        self.yVar = thermo_data

    def plot(self, filename=''):
        filename = filename if filename else "thermo_sensitivity.png"
        self.load()
        self.yVar.sort(key=lambda x: abs(x.data[-1]), reverse=True)
        if self.numSpecies:
            self.yVar = self.yVar[:self.numSpecies]
        if not self.ylabel:
            self.ylabel = 'dln(c)/d(G_i) [(kcal/mol)^-1]'
        GenericPlot.plot(self, filename=filename)

    def barplot(self, filename='', t=None):
        filename = filename if filename else "thermo_sensitivity.png"
        self.load()
        if t:
            idx = findNearest(self.xVar.data, t)
        else:
            idx = -1
        self.yVar.sort(key=lambda x: abs(x.data[idx]), reverse=True)
        if self.numSpecies:
            self.yVar = self.yVar[:self.numSpecies]

        if not self.xlabel:
            self.xlabel = 'dln(c)/d(G_i) [(kcal/mol)^-1]'
        GenericPlot.barplot(self, filename=filename, idx=idx)

    def uncertaintyPlot(self, totalVariance, t=None, filename=''):
        """
        Plot the top uncertainty contributions resulting from uncertainties in the
        thermo parameters.  The totalVariance must be specified.  Optionally,
        the reaction time `t` in seconds can be specified for plotting the uncertainties.
        The number of thermo uncertainties to plot is determined by self.numSpecies
        """

        filename = filename if filename else "thermo_uncertainty.png"
        self.load()
        if t:
            idx = findNearest(self.xVar.data, t)
        else:
            idx = -1

        thermo_uncertainty = []

        total_uncertainty = totalVariance

        for thermoSens in self.yVar:
            if isinstance(thermoSens, np.ndarray):
                # The parameter uncertainties are an array which should have the same length as the sensitivity data
                uncertainty_data = thermoSens * thermoSens.uncertainty
                uncertainty_contribution = uncertainty_data[idx] ** 2
            else:
                # The parameter uncertainty is a scalar
                uncertainty_contribution = (thermoSens.data[idx] * thermoSens.uncertainty) ** 2
            thermo_uncertainty.append([thermoSens.label, thermoSens.species, uncertainty_contribution])

        # Normalize and create new list of GenericData
        new_y_var = []
        for label, species, uncertainty in thermo_uncertainty:
            data = GenericData(label=label, species=species, data=[uncertainty / total_uncertainty * 100])
            new_y_var.append(data)

        new_y_var.sort(key=lambda x: abs(x.data[0]), reverse=True)
        new_y_var = new_y_var[:self.numSpecies]
        GenericPlot(xVar=None, yVar=new_y_var, xlabel="Uncertainty Contribution (%)").barplot(filename=filename)

        return thermo_uncertainty
