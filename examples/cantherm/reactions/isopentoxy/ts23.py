
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
    'CBS-QB3': QchemLog('ts23.out'),    
}

geometry = QchemLog('ts23.out')
frequencies = QchemLog('ts23.out')

frequencyScaleFactor = 0.983

"""pivot are the two atoms that are attached to the rotor                                                                                      
top contains the atoms that are being rotated including one of the atoms from pivots                                                           
symmetry is the symmetry number of the scan                                                                                                    
fit is fit of the scan data. It defaults to 'best', but can also be assigned as 'cosine' or 'fourier' """
rotors = [
    HinderedRotor(scanLog=QchemLog('c5h11o3scan4.out'), pivots=[1,14], top=[14,15,16,17], symmetry=3, fit='best'), #ch3, using c5h11o3scan4
    HinderedRotor(scanLog=QchemLog('c5h11o3scan4.out'), pivots=[1,10], top=[10,11,12,13], symmetry=3, fit='best'), #ch3
]

