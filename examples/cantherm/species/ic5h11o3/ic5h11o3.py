
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
    'CBS-QB3': QchemLog('c5h11o-3.out'),    
}

geometry = QchemLog('c5h11o-3.out')
frequencies = QchemLog('c5h11o-3.out')

frequencyScaleFactor = 0.983

"""pivot are the two atoms that are attached to the rotor                                                                                      
top contains the atoms that are being rotated including one of the atoms from pivots                                                           
symmetry is the symmetry number of the scan                                                                                                    
fit is fit of the scan data. It defaults to 'best', but can also be assigned as 'cosine' or 'fourier' """
rotors = [
    HinderedRotor(scanLog=QchemLog('c5h11o3scan1.out'), pivots=[1,3], top=[3, 4, 5, 6, 7, 8, 9], symmetry=1, fit='best'), #middle rot
    HinderedRotor(scanLog=QchemLog('c5h11o3scan2.out'), pivots=[3,5], top=[5,6,7,8,9], symmetry=1, fit='best'), # ch2oh 
    HinderedRotor(scanLog=QchemLog('c5h11o3scan3.out'), pivots=[1,14], top=[14,15,16,17], symmetry=3, fit='best'), #ch3
    HinderedRotor(scanLog=QchemLog('c5h11o3scan4.out'), pivots=[1,10], top=[10,11,12,13], symmetry=3, fit='best'), #ch3
    HinderedRotor(scanLog=QchemLog('c5h11o3scan5.out'), pivots=[5,8], top=[8,9], symmetry=1, fit='best'), #oh
]

