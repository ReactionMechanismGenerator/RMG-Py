"""
A list of the names of the species whose observables will be use to reduce the model.

A species observable may be the mole fraction at the end time of the batch reactor simulation.

The species name is the name of the species that can be found in the chemkin model or species 
dictionary, MINUS the species index in parentheses.

E.g.: ethane(1) --> ethane
"""

targets = ['ethane']

"""
A value between 0 and 1 that indicates how much percent (1 = 100%) the species observables
may deviate from their original values.

0.05 = 5% relative deviation
"""

tolerance = .05
