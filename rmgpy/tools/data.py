
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