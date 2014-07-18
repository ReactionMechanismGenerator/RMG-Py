
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
    'CBS-QB3': QchemLog('iC5H11O-2_betasc_CH2CHOH2.out'),    
}

geometry = QchemLog('iC5H11O-2_betasc_CH2CHOH2.out')
frequencies = QchemLog('iC5H11O-2_betasc_CH2CHOH2.out')

frequencyScaleFactor = 0.983

"""pivot are the two atoms that are attached to the rotor                                                                                      
top contains the atoms that are being rotated including one of the atoms from pivots                                                           
symmetry is the symmetry number of the scan                                                                                                    
fit is fit of the scan data. It defaults to 'best', but can also be assigned as 'cosine' or 'fourier' """
rotors = [
    HinderedRotor(scanLog=QchemLog('c5h11o2scan4.out'), pivots=[2,3], top=[3,11,12,13], symmetry=3, fit='best'), #ch3, using c5h11o2scan4.out
    HinderedRotor(scanLog=QchemLog('c5h11o2scan4.out'), pivots=[1,2], top=[1,8,9,10], symmetry=3, fit='best'), #ch3
    HinderedRotor(scanLog=QchemLog('c5h11o2scan5.out'), pivots=[5,6], top=[6,17], symmetry=1, fit='best'), # oh
]

