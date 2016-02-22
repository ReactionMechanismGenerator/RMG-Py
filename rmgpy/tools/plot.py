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
