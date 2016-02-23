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
    import numpy
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
    import numpy
    idx = (numpy.abs(array-value)).argmin()
    return idx

class GenericData(object):
    """
    A generic data class for the purpose of plotting.
    label is the original label for the data
    data is the numpy array containing the data
    The other parameters are additional common flags for data found in RMG
    """
    def __init__(self, label='', data=None, species=None, reaction=None, units=None, index=None):
        self.label = label
        self.data = data
        self.species = species
        self.reaction = reaction
        self.units = units
        self.index = index

class GenericPlot(object):
    """
    A generic plotting class that can be extended to plot other things.
    """
    def __init__(self, xVar=None, yVar=None, title='', xlabel='', ylabel=''):
        self.xVar = xVar
        self.yVar = yVar
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        
    def plot(self, filename=''):
        """
        Execute the actual plotting
        """
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        
        mpl.rc('font',family='monospace')
        fig=plt.figure()
        
        ax = fig.add_subplot(111)
        xVar = self.xVar
        yVar = self.yVar
        # Convert yVar to a list if it wasn't one already
        if isinstance(yVar, GenericData):
            yVar = [yVar]
        
            
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
            width = 1.1 + .011*maxStringLength
            legend = ax.legend(handles,labels,loc='upper center', numpoints=1, bbox_to_anchor=(width,1)) #bbox_to_anchor=(1.01,.9)
            fig.savefig(filename, bbox_extra_artists=(legend,), bbox_inches='tight')
        else:
            fig.savefig(filename)
            
    def barplot(self, filename='', idx=None):
        """
        Plot a generic barplot using just the yVars.
        idx is the index of the each y-variable to be plotted. if not given, the last value will be used
        """
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        import numpy
        
        mpl.rc('font',family='monospace')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        position = numpy.arange(len(self.yVar),0,-1)+0.5
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
            
        fig.savefig(filename)
    
    def comparePlot(self, otherGenericPlot, filename=''):
        """
        Plot a comparison data plot of this data vs a second GenericPlot class
        """
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        
        mpl.rc('font',family='monospace')
        fig=plt.figure()
        
        ax = fig.add_subplot(111)
        
        styles = ['-',':']
        # Plot the sets of data
        for i, plot in enumerate([self, otherGenericPlot]):
            # Reset the color cycle per plot to get matching colors in each set
            plt.gca().set_prop_cycle(None)
            
            xVar = plot.xVar
            yVar = plot.yVar
            # Convert yVar to a list if it wasn't one already
            if isinstance(yVar, GenericData):
                yVar = [yVar]
            
                
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
                    ax.plot(xVar.data, y.data, styles[i], label=y.label)
            
            # Plot the second set of data
            
        # Use the labels from this data object
        
        if self.xlabel:
            plt.xlabel(self.xlabel)
        elif self.xVar.label:
            xlabel = self.xVar.label
            if self.xVar.units: xlabel += ' ({0})'.format(self.xVar.units)
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
            width = 1.2+ .011*maxStringLength*2
            legend = ax.legend(handles,labels,loc='upper center', numpoints=1, bbox_to_anchor=(width,1), ncol=2) #bbox_to_anchor=(1.01,.9)
            fig.savefig(filename, bbox_extra_artists=(legend,), bbox_inches='tight')
        else:
            fig.savefig(filename)
    
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
    def __init__(self, xVar=None, yVar=None, title='', xlabel='', ylabel='', csvFile='', numSpecies=None, species={}):
        GenericPlot.__init__(self, xVar=xVar, yVar=yVar, title=title, xlabel=xlabel, ylabel=ylabel)
        self.csvFile = csvFile
        self.numSpecies = numSpecies
        self.species = species
        
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
    
    
    def comparePlot(self, otherSimulationPlot, filename=''):
        
        filename = filename if filename else 'simulation_compare.png'
        self.load()
        otherSimulationPlot.load()
        GenericPlot.comparePlot(self, otherSimulationPlot, filename)
        
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
    def __init__(self, xVar=None, yVar=None, title='', xlabel='', ylabel='', csvFile='', numReactions=None, reactions={}):
        GenericPlot.__init__(self, xVar=xVar, yVar=yVar, title=title, xlabel=xlabel, ylabel=ylabel)
        self.csvFile = csvFile
        self.numReactions = numReactions
        self.reactions = reactions
        
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
        self.ylabel = 'dlnc/dlnk'
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
        
        self.xlabel = 'dlnc/dlnk_i'
        GenericPlot.barplot(self, filename=filename, idx=idx)
        
class ThermoSensitivityPlot(GenericPlot):
    """
    A class for plotting the top sensitivities to a thermo DeltaG value of species within the model.
    The value used is the sensitivity at the final time point.
    `numSpecies` indicates the number of species to plot. 
    
    `species` is a dictionary corresponding to specific species thermo sensitivities to be plotted
    
    barplot() will instead plot a horizontal bar plot of the sensitivities at a given
    time step.  If time step is not given, the end step will automatically be chosen
    """
    def __init__(self, xVar=None, yVar=None, title='', xlabel='', ylabel='', csvFile='', numSpecies=None, species={}):
        GenericPlot.__init__(self, xVar=xVar, yVar=yVar, title=title, xlabel=xlabel, ylabel=ylabel)
        self.csvFile = csvFile
        self.numSpecies = numSpecies
        self.species = species
    
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
        self.ylabel = 'dlnc/dlnG'
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
        self.xlabel = 'dlnc/dlnG'
        GenericPlot.barplot(self, filename=filename, idx=idx)