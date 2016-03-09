import numpy

class GenericData(object):
    """
    A generic data class for the purpose of plotting.
    ======================= ==============================================================================================
    Attribute               Description
    ======================= ==============================================================================================
    `label`                 A string label describing the data, can be used in a plot legend or in an axis label
    `data`                  A numpy array of the data
    `species`               Contains species associated with the data, often used with a Species object
    `reaction`              Contains reaction associated with the data, often used with a Reaction object
    `units`                 Contains a string describing the units associated with the data
    `index`                 An integer containing the index associated with the data
    ======================= ==============================================================================================
    """
    def __init__(self, label='', data=None, species=None, reaction=None, units=None, index=None):
        
        self.label = str(label) if label else None
        
        if isinstance(data, list):
                self.data = numpy.array(data)
        elif isinstance(data, numpy.ndarray):
                self.data = data
        else:
            raise Exception('Data for GenericData object must be initialized as a list or numpy.array of values.')
            
        self.species = species
        self.reaction = reaction
        self.units = str(units) if units else None
        self.index = int(index) if index else None