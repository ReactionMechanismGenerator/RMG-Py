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

import matplotlib as mpl
# Force matplotlib to not use any Xwindows backend.
# This must be called before pylab, matplotlib.pyplot, or matplotlib.backends is imported
mpl.use('Agg')
import matplotlib.pyplot as plt
from rmgpy.tools.data import GenericData
import numpy
import csv
import math
        
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
    indexPattern = re.compile(r'^\S+\(\d+\)$')
    unitsPattern = re.compile(r'\s\(.+\)$')
    rxnSensPattern = re.compile('^dln\[\S+\]\/dln\[k\d+\]:\s\S+$')
    thermoSensPattern = re.compile('^dln\[\S+\]\/dG\[\S+\]$')
    
    timeData = []; data = {}
    f = csv.reader(open(csvFile, 'r'))
    
    columns = zip(*f)
    time = GenericData(label = columns[0][0],
                       data = numpy.array(columns[0][1:],dtype=numpy.float64),
                      )
    
    # Parse the units from the Time header
    if unitsPattern.search(time.label):
        label, sep, units = time.label[:-1].rpartition('(')
        time.label = label
        time.units = units
            
            
    dataList = []
    for col in columns[1:]:
        header = col[0]
        values = numpy.array(col[1:],dtype=numpy.float64)
        data = GenericData(label=header,data=values)
        
        # Parse the index or the label from the header
        if indexPattern.search(data.label):
            species, sep, index = data.label[:-1].rpartition('(')
            # Save the species attribute if an index was found
            data.species = species
            data.index = int(index)
        elif unitsPattern.search(data.label):
            label, sep, units = data.label[:-1].rpartition('(')
            data.label = label
            data.units = units
        elif rxnSensPattern.search(data.label):
            rxn = data.label.split()[1]
            index = data.label.split()[0][:-2].rpartition('dln[k')[2]
            data.reaction = rxn
            data.index = index
        elif thermoSensPattern.search(data.label):
            species = data.label[:-1].rpartition('dG[')[2]
            data.species = species
            if indexPattern.search(species):
                data.index = species[:-1].rpartition('(')[2]
            
        dataList.append(data)
        
        
    return time, dataList

def findNearest(array, value):
    """
    Returns the index of the closest value in a sorted array
    """
    idx = (numpy.abs(array-value)).argmin()
    return idx

def linearlyInterpolatePoint(xArray, yArray, xValue):
    """
    Returns the interpolated yValue for given xValue using data from the two sorted arrays:
    """
    #Find the next largest point in xArray that is still smaller than xValue:
    lowerIndex=None
    for index, x in enumerate(xArray):
        if x>xValue:
            break
        lowerIndex=index

    #If xValue is outside the domain of xArray, we use either the min or max points for dydx
    if lowerIndex is None:
        lowerIndex=0
    elif lowerIndex==len(xArray)-1:
        lowerIndex=lowerIndex-1
    higherIndex=lowerIndex+1

    dydx=(yArray[higherIndex]-yArray[lowerIndex])/(xArray[higherIndex]-xArray[lowerIndex])

    if xValue < xArray[lowerIndex]:
        yValue=yArray[lowerIndex]-dydx*(xValue-xArray[lowerIndex])
    else:
        yValue=yArray[lowerIndex]+dydx*(xValue-xArray[lowerIndex])
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
        mpl.rc('font',family='sans-serif')
        fig=plt.figure()
        
        ax = fig.add_subplot(111)

        xVar = self.xVar
        yVar = self.yVar
            
        if len(yVar) == 1:
            y = yVar[0]
            ax.plot(xVar.data, y.data)
            
            # Create a ylabel for the label of the y variable
            if not self.ylabel and y.label:
                ylabel = y.label
                if y.units: ylabel += ' ({0})'.format(y.units)
                plt.ylabel(ylabel)
        else:
            for y in yVar:
                ax.plot(xVar.data, y.data, '-', label=y.label)
        
        if self.xlabel:
            plt.xlabel(self.xlabel)
        elif xVar.label:
            xlabel = xVar.label
            if xVar.units: xlabel += ' ({0})'.format(xVar.units)
            plt.xlabel(xlabel)
            
        if self.ylabel:
            plt.ylabel(self.ylabel)
            
        if self.title:
            plt.title(self.title)
            
        ax.grid('on')
        handles, labels = ax.get_legend_handles_labels()
        if labels:
            # Create a legend outside the plot and adjust width based off of longest legend label
            maxStringLength = max([len(label) for label in labels])
            width = 1.05 + .011*maxStringLength
            legend = ax.legend(handles,labels,loc='upper center', numpoints=1, bbox_to_anchor=(width,1)) #bbox_to_anchor=(1.01,.9)
            fig.savefig(filename, bbox_extra_artists=(legend,), bbox_inches='tight')
        else:
            fig.savefig(filename, bbox_inches='tight')
            
    def barplot(self, filename='', idx=None):
        """
        Plot a generic barplot using just the yVars.
        idx is the index of the each y-variable to be plotted. if not given, the last value will be used
        """
        mpl.rc('font',family='sans-serif')
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        position = numpy.arange(len(self.yVar),0,-1)
        # Reverse in order to go front top to bottom
        if not idx:
            idx = -1
        ax.barh(position, numpy.array([y.data[idx] for y in self.yVar]), align='center', alpha=0.5)
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
    
    def comparePlot(self, otherGenericPlots, filename='', title='', xlabel='', ylabel='', styles=None):
        """
        Plot a comparison data plot of this data vs a other Generic Plot classes

        otherGenericPlots is a single GenericPlot or list of GenericPlots to plot against

        styles is a list of strings corresponding to matplotlib shorthand line/marker styles for each corresponding
            other Generic Plot.
            -If the number of plots is greater than 4, then it is mandatory to give a custom styles
        """
        
        mpl.rc('font',family='sans-serif')
        #mpl.rc('text', usetex=True)
        fig=plt.figure()
        
        ax = fig.add_subplot(111)

        if styles is None:
            styles = ['-','--', ':', '-.']
        if not type(otherGenericPlots) is list:
            otherGenericPlots=[otherGenericPlots]
        # Plot the sets of data
        for i, plot in enumerate([self]+otherGenericPlots):
            # Reset the color cycle per plot to get matching colors in each set
            plt.gca().set_prop_cycle(None)
            
            xVar = plot.xVar
            yVar = plot.yVar
            # Convert yVar to a list if it wasn't one already
            if isinstance(yVar, GenericData):
                yVar = [yVar]
            
                
            if len(yVar) == 1:
                y = yVar[0]
                ax.plot(xVar.data, y.data, styles[i], label=y.label)
                # Save a ylabel based on the y variable's label if length of yVar contains only 1 variable
                if not self.ylabel and y.label:
                    self.ylabel = y.label
                    if y.units: self.ylabel += ' ({0})'.format(y.units)
            else:
                for y in yVar:
                    ax.plot(xVar.data, y.data, styles[i], label=y.label)
            
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
            
        ax.grid('on')
        handles, labels = ax.get_legend_handles_labels()
        if labels:
            # Create a legend outside the plot and adjust width based off of longest legend label
            maxStringLength = max([len(label) for label in labels])
            width = 1.2+ .011*maxStringLength*2
            legend = ax.legend(handles,labels,loc='upper center', numpoints=1, bbox_to_anchor=(width,1), ncol=2) #bbox_to_anchor=(1.01,.9)
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
            time, dataList = parseCSVData(self.csvFile)
        else:
            time = self.xVar
            dataList = self.yVar
            
        speciesData = []
        if self.species:
            # A specific set of species was specified to be plotted
            for speciesLabel, chemkinLabel in self.species.iteritems():
                for data in dataList:
                    if chemkinLabel == data.label:
                        # replace the data label with the desired species label
                        data.label = speciesLabel
                        speciesData.append(data)
                        break 
        
        else:
            for data in dataList:
                # Only plot if RMG detects that the data corresponds with a species
                # This will not include bath gases
                if data.species:
                    speciesData.append(data)
            
        self.xVar = time
        self.yVar = speciesData
            
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
    def __init__(self, xVar=None, yVar=None, title='', xlabel='', ylabel='', csvFile='', numReactions=None, reactions=None):
        GenericPlot.__init__(self, xVar=xVar, yVar=yVar, title=title, xlabel=xlabel, ylabel=ylabel)
        self.csvFile = csvFile
        self.numReactions = numReactions
        self.reactions = reactions if reactions else {}
        
    def load(self):
        if self.xVar == None and self.yVar == None:
            time, dataList = parseCSVData(self.csvFile)
        else:
            time = self.xVar
            dataList = self.yVar
        reactionData = []
        if self.reactions:
            # A specific set of reaction sensitivities was specified to be plotted
            for reactionLabel, chemkinLabel in self.reactions.iteritems():
                for data in dataList:
                    if chemkinLabel == data.reaction:
                        # replace the data label with the desired species label
                        data.label = reactionLabel
                        reactionData.append(data)
                        break 
        else:
            for data in dataList:
                if data.reaction:
                    reactionData.append(data)
        self.xVar = time
        self.yVar = reactionData
            
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

        reactionUncertainty = []

        totalUncertainty = totalVariance

        for reactionSens in self.yVar:
            if isinstance(reactionSens,numpy.ndarray):
                # The parameter uncertainties are an array which should have the same length as the sensitivity data
                uncertaintyData = reactionSens*reactionSens.uncertainty
                uncertaintyContribution = uncertaintyData[idx]**2
            else:
                uncertaintyContribution = (reactionSens.data[idx]*reactionSens.uncertainty)**2

            reactionUncertainty.append([reactionSens.label, reactionSens.reaction, uncertaintyContribution])

        # Normalize and create new list of GenericData
        newYVar = []
        for label, reaction, uncertainty in reactionUncertainty:
            data = GenericData(label=label, reaction=reaction, data = [uncertainty/totalUncertainty*100])
            newYVar.append(data)

        newYVar.sort(key=lambda x: abs(x.data[0]), reverse = True)
        newYVar = newYVar[:self.numReactions]

        GenericPlot(xVar=None, yVar=newYVar, xlabel ="Uncertainty Contribution (%)").barplot(filename=filename)
        
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
            time, dataList = parseCSVData(self.csvFile)
        else:
            time = self.xVar
            dataList = self.yVar
        thermoData = []
        if self.species:
            # A specific set of species sensitivities was specified to be plotted
            for speciesLabel, chemkinLabel in self.species.iteritems():
                for data in dataList:
                    if chemkinLabel == data.species:
                        # replace the data label with the desired species label
                        data.label = speciesLabel
                        thermoData.append(data)
                        break 
        else:
            for data in dataList:
                if data.species:
                    thermoData.append(data)
                    
        self.xVar = time
        self.yVar = thermoData
            
    def plot(self, filename=''):
        filename = filename if filename else "thermo_sensitivity.png"
        self.load()
        self.yVar.sort(key=lambda x: abs(x.data[-1]), reverse = True)
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
        self.yVar.sort(key=lambda x: abs(x.data[idx]), reverse = True)
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

        thermoUncertainty = []

        totalUncertainty = totalVariance

        for thermoSens in self.yVar:
            if isinstance(thermoSens,numpy.ndarray):
                # The parameter uncertainties are an array which should have the same length as the sensitivity data
                uncertaintyData = thermoSens*thermoSens.uncertainty
                uncertaintyContribution = uncertaintyData[idx]**2
            else:
                # The parameter uncertainty is a scalar
                uncertaintyContribution = (thermoSens.data[idx]*thermoSens.uncertainty)**2
            thermoUncertainty.append([thermoSens.label, thermoSens.species, uncertaintyContribution])


        # Normalize and create new list of GenericData
        newYVar = []
        for label, species, uncertainty in thermoUncertainty:
            data = GenericData(label=label, species = species, data = [uncertainty/totalUncertainty*100])
            newYVar.append(data)


        newYVar.sort(key=lambda x: abs(x.data[0]), reverse = True)
        newYVar = newYVar[:self.numSpecies]
        GenericPlot(xVar=None, yVar=newYVar, xlabel ="Uncertainty Contribution (%)").barplot(filename=filename)

class ComparisonBundle:
    """
    A class for storing multiple :class: GenericData for the purpose of comparison
    ======================= ==============================================================================================
    Attribute               Description
    ======================= ==============================================================================================
    `title`                 A string label describing the data, (recommended to use dependent variable and use source
                                differentiation for label of GenericData objects in yData)
    `xDataList`             A list of :class: GenericData for the independent variable (list index corresponds to yData)
    'yDataList'             A list of :class: GenericData the independent variable (list index corresponds to xData)
    `species`               Contains species associated with the data, often used with a Species object
    `reaction`              Contains reaction associated with the data, often used with a Reaction object
    `xUnits`                Contains a string describing the units associated with the data of the independent variable
    'yUnits'                Contains a string describing the units associated with the data of the independent variable
    `index`                 An integer containing the index associated with the data
    ======================= ==============================================================================================
    """
    def __init__(self, title='', xDataList=None, yDataList=None, species=None, reaction=None, xUnits=None, yUnits=None, index=None):
        self.title=title
        if xDataList:
            self.xDataList=xDataList
        else:
            self.xDataList=[]
        if yDataList:
            self.yDataList=yDataList
        else:
            self.yDataList=[]
        self.species=species
        self.reaction=reaction
        self.xUnits=xUnits
        self.yUnits=yUnits
        self.index=index

        #Check that there is a xData for every yData
        assert len(self.xDataList) == len(self.yDataList), "The length of xDataList and yDataList are not the same."

        #assign indicies to the xData and yData
        for index, (x,y) in enumerate(zip(self.xDataList, self.yDataList)):
            x.index=index
            y.index=index

        #Check that species, reaction, and unit are consistent across all data
        self.checkAndMakeConsistent()

    def __str__(self):
        """
        Return a string representation of this test case, using its title'.
        """
        return 'Comparison Bundle: {0}'.format(self.title)

    def checkAndMakeConsistent(self):
        """
        Checks that the species, reaction, and units are consistent across all data, raising an assertion error if
         anything with the exception of an omission in the form of None.

        If species, reaction, yUnits, or xUnits are None on the ComparisonBundle (self), then this function will
        attempt to pick it out from GenericData objects in xDataList and yDataList

        Conversely, if species, reaction, or units are None in any of the GenericData objects in xDataList or
        yDataList, this function will set those attributes equal to the appropriate analog from the ComparisonBundle.
        """

        #Check that there is a xData for every yData
        assert len(self.xDataList) == len(self.yDataList), "The length of xDataList and yDataList are not the same."

        #If there is no data, there is nothing to check
        if not self.xDataList>0:
            return

        #Matching of attributes for the "head" (self) and the genericDatas
        matching={'species': (self.yDataList, 'species', "yData"),
               'reaction': (self.yDataList, 'reaction', "xData"),
               'xUnits': (self.xDataList, 'units', 'xData'),
               'yUnits': (self.yDataList, 'units', 'yData'),

        }

        #If the head attributes are not set, try to take it from one of the arrays in the list
        for headAttr, objAttr in matching.iteritems():
            dataList=objAttr[0]
            slaveAttr=objAttr[1]
            if getattr(self, headAttr) is None:
                for data in dataList:
                    if not getattr(data, slaveAttr) is None:
                        setattr(self, headAttr, getattr(data, slaveAttr))

        #Now check that all attributes are consistent
        for headAttr, objAttr in matching.iteritems():
            dataList=objAttr[0]
            slaveAttr=objAttr[1]
            for data in dataList:
                if getattr(data, slaveAttr) is None:
                    setattr(data, slaveAttr, getattr(self, headAttr))

                if getattr(self,headAttr) is None:
                    assert getattr(data, slaveAttr)==getattr(self, headAttr), \
                        "The GenericData for {0} index {1} has inconsistent {2} with the '{3}' or " \
                        "other GenericData objects in the ComparisonBundle".format(objAttr[2], data.index, slaveAttr,
                                                                                   str(self))
                elif headAttr== "species" or headAttr=="reaction":
                    assert getattr(self, headAttr).isIsomorphic(getattr(data, slaveAttr)), \
                        "The GenericData for {0} index {1} has inconsistent {2} with the '{3}' or " \
                        "other GenericData objects in the ComparisonBundle".format(objAttr[2], data.index, slaveAttr,
                                                                                   str(self))
                else:
                    assert getattr(data, slaveAttr)==getattr(self, headAttr), \
                        "The GenericData for {0} index {1} has inconsistent {2} with the '{3}' or " \
                        "other GenericData objects in the ComparisonBundle".format(objAttr[2], data.index, slaveAttr,
                                                                                   str(self))
    def addDataSet(self, xData, yData):
        """
        Adds a new set of data to the ComparisonBundle where xData and yData are :class: GenericData objects with either
        matching species, reaction, and units, or with all the former units unintialized with None
        """
        newIndex=len(self.xDataList)
        xData.index=newIndex
        self.xDataList.append(xData)

        yData.index=newIndex
        self.yDataList.append(yData)

        #Recheck that the newly added data set conforms to the correct attributes
        self.checkAndMakeConsistent()

    def removeDataSet(self, yLabel=None, index=None):
        """
        Removes a data set (a pair of :class: GenericData objects, one from xDataList and one from yDataList). The data
        set is recognized by the label of the yData OR the index of data set which is the same as the list index in
        xDataList or yDataList.
        """
        if index is None and yLabel is None:
            raise Exception("No label or index was given to identify the Data set to be removed")

        #Find index if label is provided
        if yLabel:
            for newIndex, yData in enumerate(self.yDataList):
                if yData.label==yLabel:
                    index=newIndex
                    break
            else: raise Exception("The inputted label {0} did not match the label for any of the data sets " \
                                  "(checks the y variable)".format(yLabel))

        #Remove the data set
        del self.xDataList[index]
        del self.yDataList[index]

        #reindex all remaining data sets
        for xData, yData in zip(self.xDataList[index:], self.yDataList[index:]):
            xData.index+=-1
            yData.index+=-1

    def makeLog(self, base=10.0):
        """
        Converts the yData arrays to logarithms. The variable 'base' is a float defining the base of the logarithm.

        Any data points that are not compatible (negative value or 0), will be discarded along with the corresponding
        point in xData
        """

        for xData, yData in zip(self.xDataList, self.yDataList):
            newXArray=[]
            newYArray=[]
            for xPoint, yPoint in zip(xData.data, yData.data):
                if yPoint >0:
                    newXArray.append(xPoint)
                    newYArray.append(math.log(yPoint,base))
            #Convert to numpy array and set data
            newXArray=numpy.array(newXArray)
            newYArray=numpy.array(newYArray)
            xData.data=newXArray
            yData.data=newYArray

            #Change units on yData
            yData.units="log" +str(base)+ " " + yData.units

        #Change units on ComparisonBundle
        self.yUnits="log" +str(base)+ " " + self.yUnits

    def sortByX(self):
        """
        Sort all the datapoints in each GenericData in xDataList in ascending order and perform a corresponding sort
        on the datapoints in each GenericData in yDataList
        """

        for xData, yData in zip(self.xDataList, self.yDataList):
            #decorated sort
            xTuple, yTuple=zip(*sorted(zip(xData.data, yData.data)))
            #set data
            xData.data=numpy.array(xTuple)
            yData.data=numpy.array(yTuple)

    def printToCsv(self, outputPath):
        """
        Args:
            outputPath: Pathway to where csv will be saved

        Returns: Saves a csv file of human readable data for the ComparisonBundle

        """
        #Sort to make human readable
        self.sortByX()
        with open(outputPath, 'wb') as output:
            spamwriter=csv.writer(output)
            #write out headings and units
            labelLine=[]
            unitLine=[]
            for xData, yData in zip(self.xDataList, self.yDataList):
                labelLine.append(xData.label)
                labelLine.append(yData.label)
                unitLine.append(xData.units)
                unitLine.append(yData.units)

            spamwriter.writerow(labelLine)
            spamwriter.writerow(unitLine)

            dataMatrix=[]
            largestLength=max([len(xData.data) for xData in self.xDataList])
            for xData, yData in zip(self.xDataList, self.yDataList):
                #Concatenate extra blank entries if length is not longest
                if len(xData.data)<largestLength:
                    tail=numpy.array(['']*(largestLength-len(xData.data)))
                    newXData=numpy.concatenate((xData.data,tail))
                    newYData=numpy.concatenate((yData.data,tail))
                    dataMatrix.append(newXData)
                    dataMatrix.append(newYData)
                else:
                    dataMatrix.append(xData.data)
                    dataMatrix.append(yData.data)

            #transpose dataMatrix
            numpy.array(dataMatrix)
            newDataMatrix=numpy.transpose(dataMatrix)

            #write out data points
            for line in newDataMatrix:
                spamwriter.writerow(line)

    def plot(self, filename):
        """
        Function to plot a Comparison Bundle

        filename is str path to where the plot.png will be saved
        """
        #Take title, xlabel and ylabel from ComparisonBundle
        title= self.title + " Comparison"
        xlabel= "{0} ({1})".format(self.xDataList[0].label, self.xUnits)
        ylabel= "{0} ({1})".format(self.title, self.yUnits)

        #Make Generic Plot objects out of xData and yData
        plotList=[]

        for xData, yData in zip(self.xDataList, self.yDataList):
            #Careful that the ylabel here yData.label (so each line has a separate legend entry) instead of ylabel
            plotList.append(GenericPlot(xData, yData, title, xlabel, yData.label))

        plotList[0].comparePlot(plotList[1:], filename, ylabel=ylabel, styles=['o', '-', '--'])

