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