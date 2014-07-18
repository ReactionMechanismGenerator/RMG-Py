
#!/usr/bin/env python                                                                                                                          
# encoding: utf-8

atoms = {
    'C': 5,
    'H': 11,
    'O': 1
}

bonds = {}
linear = False

externalSymmetry = 1

spinMultiplicity = 2

opticalIsomers = 1

energy = {
    'CBS-QB3': QchemLog('ts24.out'),    
}

geometry = QchemLog('ts24.out')
frequencies = QchemLog('ts24.out')

frequencyScaleFactor = 0.983

"""pivot are the two atoms that are attached to the rotor                                                                                      
top contains the atoms that are being rotated including one of the atoms from pivots                                                           
symmetry is the symmetry number of the scan                                                                                                    
fit is fit of the scan data. It defaults to 'best', but can also be assigned as 'cosine' or 'fourier' """
rotors = [
    HinderedRotor(scanLog=QchemLog('c5h11o4scan3.out'), pivots=[1,6], top=[6,7,8,9], symmetry=3, fit='best'), #ch3
    HinderedRotor(scanLog=QchemLog('c5h11o4scan4.out'), pivots=[1,2], top=[2,3,4,5], symmetry=3, fit='best'), #ch3
]

