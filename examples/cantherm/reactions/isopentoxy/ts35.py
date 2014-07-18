
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
    'CBS-QB3': QchemLog('ts35.out'),    
}

geometry = QchemLog('ts35.out')
frequencies = QchemLog('ts35.out')

frequencyScaleFactor = 0.983

"""pivot are the two atoms that are attached to the rotor                                                                                      
top contains the atoms that are being rotated including one of the atoms from pivots                                                           
symmetry is the symmetry number of the scan                                                                                                    
fit is fit of the scan data. It defaults to 'best', but can also be assigned as 'cosine' or 'fourier' """
rotors = [
    HinderedRotor(scanLog=QchemLog('c5h11o5scan4.out'), pivots=[3,4], top=[4,10], symmetry=1, fit='best'), #oh
    HinderedRotor(scanLog=QchemLog('c5h11o5scan5.out'), pivots=[2,3], top=[3,4,9,10,17], symmetry=1, fit='best'), #ch2oh
    HinderedRotor(scanLog=QchemLog('c5h11o5scan2.out'), pivots=[1,6], top=[6,14,15,16], symmetry=3, fit='best'), #ch3
]

