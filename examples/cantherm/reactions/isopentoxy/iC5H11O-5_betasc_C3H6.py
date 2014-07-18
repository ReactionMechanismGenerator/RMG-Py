
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
    'CBS-QB3': QchemLog('iC5H11O-5_betasc_C3H6m08.out'),    
}

geometry = QchemLog('iC5H11O-5_betasc_C3H6m08.out')
frequencies = QchemLog('iC5H11O-5_betasc_C3H6m08.out')

frequencyScaleFactor = 0.983

"""pivot are the two atoms that are attached to the rotor                                                                                      
top contains the atoms that are being rotated including one of the atoms from pivots                                                           
symmetry is the symmetry number of the scan                                                                                                    
fit is fit of the scan data. It defaults to 'best', but can also be assigned as 'cosine' or 'fourier' """
rotors = [
    HinderedRotor(scanLog=QchemLog('c5h11o5scan2.out'), pivots=[2,3], top=[3,7,8,9], symmetry=3, fit='best'), #ch3
    HinderedRotor(scanLog=QchemLog('c5h11o5scan4.out'), pivots=[12,11], top=[12,17], symmetry=1, fit='best'), #oh
#    HinderedRotor(scanLog=QchemLog('c5h11o5scan5.out'), pivots=[10,13], top=[13,14,15,16,17], symmetry=1, fit='best'), #ch2oh
]

