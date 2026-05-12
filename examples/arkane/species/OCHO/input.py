title = 'OCHO'

description = \
"""
An example for defining a species using a direct (numeric) E0 input (see in OCHO.py)
"""


species('OCHO', 'OCHO.py',
        structure = SMILES('[O]C=O'),
        collisionModel = TransportData(sigma=(3.59,'angstrom'), epsilon=(4140.61,'J/mol')),
        energyTransferModel = SingleExponentialDown(alpha0=(80,'cm^-1'), T0=(300,'K'), n=0.85),
)

statmech('OCHO')

