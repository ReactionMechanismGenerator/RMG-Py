#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
This file provides atomic energy and BAC parameters for several model
chemistries.
"""

# Atom energy corrections to reach gas-phase reference state
# Experimental enthalpy of formation at 0 K, 1 bar for gas phase
# See Gaussian thermo whitepaper at http://gaussian.com/thermo/
# Note: These values are relatively old and some improvement may be possible by using newer values
# (particularly for carbon).
# However, care should be taken to ensure that they are compatible with the BAC values (if BACs are used)
# He, Ne, K, Ca, Ti, Cu, Zn, Ge, Br, Kr, Rb, Ag, Cd, Sn, I, Xe, Cs, Hg, and Pb are taken from CODATA
# Codata: Cox, J. D., Wagman, D. D., and Medvedev, V. A., CODATA Key Values for Thermodynamics, Hemisphere
# Publishing Corp., New York, 1989. (http://www.science.uwaterloo.ca/~cchieh/cact/tools/thermodata.html)
atom_hf = {'H': 51.633, 'He': -1.481,
           'Li': 37.69, 'Be': 76.48, 'B': 136.2, 'C': 170.0284, 'N': 112.471, 'O': 58.997, 'F': 18.464, 'Ne': -1.481,
           'Na': 25.69, 'Mg': 34.87, 'Al': 78.23, 'Si': 106.6, 'P': 75.42, 'S': 65.66, 'Cl': 28.59,
           'K': 36.841, 'Ca': 41.014, 'Ti': 111.2, 'Cu': 79.16, 'Zn': 29.685, 'Ge': 87.1, 'Br': 28.182,
           'Kr': -1.481,
           'Rb': 17.86, 'Ag': 66.61, 'Cd': 25.240, 'Sn': 70.50, 'I': 24.04, 'Xe': -1.481,
           'Cs': 16.80, 'Hg': 13.19, 'Pb': 15.17}

# Thermal contribution to enthalpy for the atoms reported by Gaussian thermo whitepaper
# This will be subtracted from the corresponding value in atom_hf to produce an enthalpy used in calculating
# the enthalpy of formation at 298 K
atom_thermal = {'H': 1.01, 'He': 1.481,
                'Li': 1.1, 'Be': 0.46, 'B': 0.29, 'C': 0.25, 'N': 1.04, 'O': 1.04, 'F': 1.05, 'Ne': 1.481,
                'Na': 1.54, 'Mg': 1.19, 'Al': 1.08, 'Si': 0.76, 'P': 1.28, 'S': 1.05, 'Cl': 1.1,
                'K': 1.481, 'Ca': 1.481, 'Ti': 1.802, 'Cu': 1.481, 'Zn': 1.481, 'Ge': 1.768, 'Br': 1.481,
                'Kr': 1.481,
                'Rb': 1.481, 'Ag': 1.481, 'Cd': 1.481, 'Sn': 1.485, 'I': 1.481, 'Xe': 1.481,
                'Cs': 1.481, 'Hg': 1.481, 'Pb': 1.481}

# Spin orbit correction (SOC) in Hartrees
# Values taken from ref 22 of http://dx.doi.org/10.1063/1.477794 and converted to Hartrees
# Values in milli-Hartree are also available (with fewer significant figures) from table VII of
# http://dx.doi.org/10.1063/1.473182
# Iodine SOC calculated as a weighted average of the electronic spin splittings of the lowest energy state.
# The splittings are obtained from Huber, K.P.; Herzberg, G., Molecular Spectra and Molecular Structure. IV.
# Constants of Diatomic Molecules, Van Nostrand Reinhold Co., 1979
# Halogen Atoms taken from https://cccbdb.nist.gov/elecspin.asp
SOC = {
    'H': 0.0, 'N': 0.0, 'O': -0.000355, 'C': -0.000135, 
    'S': -0.000893, 'P': 0.0, 'I': -0.011547226,
    'Cl':-0.001338, 'Br':-0.005597, 'F':-0.000614
}

# Atomic energies
# All model chemistries here should be lower-case because the user input is changed to lower-case
atom_energies = {
    # Note: If your model chemistry does not include spin orbit coupling, you should add the corrections
    # to the energies here

    # cbs-qb3 and cbs-qb3-paraskevas have the same corrections
    'cbs-qb3': {
        'H': -0.499818 + SOC['H'], 'N': -54.520543 + SOC['N'], 'O': -74.987624 + SOC['O'],
        'C': -37.785385 + SOC['C'], 'P': -340.817186 + SOC['P'], 'S': -397.657360 + SOC['S']
    },
    'cbs-qb3-paraskevas': {
        'H': -0.499818 + SOC['H'], 'N': -54.520543 + SOC['N'], 'O': -74.987624 + SOC['O'],
        'C': -37.785385 + SOC['C'], 'P': -340.817186 + SOC['P'], 'S': -397.657360 + SOC['S']
    },

    'm06-2x/cc-pvtz': {
        'H': -0.498135 + SOC['H'], 'N': -54.586780 + SOC['N'], 'O': -75.064242 + SOC['O'],
        'C': -37.842468 + SOC['C'], 'P': -341.246985 + SOC['P'], 'S': -398.101240 + SOC['S'],
        'Cl': -460.140693 + SOC['Cl'], 'F': -99.730200, 'Br': -2574.198266
    },

    'g3': {
        'H': -0.5010030, 'N': -54.564343, 'O': -75.030991, 'C': -37.827717, 'P': -341.116432, 'S': -397.961110
    },

    # 'g4': {
    #     'Br': -2573.583929 + SOC['Br'], 'C': -37.832752 + SOC['C'],'Cl': -460.013636 + SOC['Cl'],
    #     'F': -99.703566 + SOC['F'],'H': -0.500004, 'N': -54.572250, 'O': -75.044085 + SOC['O']
    # },

    'g4': {
        'Br': -2573.58537, 'C': -37.83417, 'Cl': -460.01505, 
        'F': -99.70498, 'H': -0.50142, 'N': -54.57367, 'O': -75.04550
    },

    'b3lyp-gd3bj_jun-cc-pvtz': {
        'Br': -2574.1875226196235 + SOC['Br'],'C': -37.85905670561903 + SOC['C'],'Cl': -460.1751342594461 + SOC['Cl'],
        'F': -99.76599874142578 + SOC['F'], 'O': -75.09405247620397 + SOC['O']
    },

    'm062x-gd3_jun-cc-pvtz': {
        'Br': -2574.1981614882384 + SOC['Br'], 'C': -37.84283423642995 + SOC['C'], 'Cl': -460.14102245588464 + SOC['Cl'],
        'F': -99.73219425452456 + SOC['F'], 'H': -0.49813588717679524, 'N': -54.587423613280166,
        'O': -75.06579209968129 + SOC['O']
    },

    'm06hf-gd3_jun-cc-pvtz': {
        'Br': -2574.385001963926 + SOC['Br'], 'C': -37.85704620758064 + SOC['C'],'Cl': -460.1530080643251 + SOC['Cl'],
        'F': -99.74460663490936 + SOC['F'],'O': -75.07808694358144 + SOC['O']
    },

    'mn15_jun-cc-pvtz': {
        'Br': -2574.488328990481 + SOC['Br'], 'C': -37.80974876153511 + SOC['C'], 'Cl': -460.0507670586289 + SOC['Cl'],
        'F': -99.68176734398618 + SOC['F'], 'H': -0.49850400650689464, 'N': -54.54198264839306,
        'O': -75.01779384012698 + SOC['O']
    },

    'wb97xd_jun-cc-pvtz': {
        'C': -37.84201073383711 + SOC['C'], 'O': -75.06821605326587 + SOC['O']
    },

    # * indicates that the grid size used in the [QChem] electronic
    # structure calculation utilized 75 radial points and 434 angular points
    # (i.e,, this is specified in the $rem section of the [qchem] input file as: XC_GRID 000075000434)
    'm08so/mg3s*': {
        'H': -0.5017321350 + SOC['H'], 'N': -54.5574039365 + SOC['N'],
        'O': -75.0382931348 + SOC['O'], 'C': -37.8245648740 + SOC['C'],
        'P': -341.2444299005 + SOC['P'], 'S': -398.0940312227 + SOC['S']
    },

    'klip_1': {
        'H': -0.50003976 + SOC['H'], 'N': -54.53383153 + SOC['N'], 'O': -75.00935474 + SOC['O'],
        'C': -37.79266591 + SOC['C']
    },

    # Klip QCI(tz,qz)
    'klip_2': {
        'H': -0.50003976 + SOC['H'], 'N': -54.53169400 + SOC['N'], 'O': -75.00714902 + SOC['O'],
        'C': -37.79060419 + SOC['C']
    },

    # Klip QCI(dz,tz)
    'klip_3': {
        'H': -0.50005578 + SOC['H'], 'N': -54.53128140 + SOC['N'], 'O': -75.00356581 + SOC['O'],
        'C': -37.79025175 + SOC['C']
    },

    # Klip CCSD(T)(tz,qz)
    'klip_2_cc': {
        'H': -0.50003976 + SOC['H'], 'O': -75.00681155 + SOC['O'], 'C': -37.79029443 + SOC['C']
    },

    'ccsd(t)-f12/cc-pvdz-f12_h-tz': {
        'H': -0.499946213243 + SOC['H'], 'N': -54.526406291655 + SOC['N'],
        'O': -74.995458316117 + SOC['O'], 'C': -37.788203485235 + SOC['C']
    },

    'ccsd(t)-f12/cc-pvdz-f12_h-qz': {
        'H': -0.499994558325 + SOC['H'], 'N': -54.526406291655 + SOC['N'],
        'O': -74.995458316117 + SOC['O'], 'C': -37.788203485235 + SOC['C']
    },

    # We are assuming that SOC is included in the Bond Energy Corrections
    'ccsd(t)-f12/cc-pvdz-f12': {
        'H': -0.499811124128, 'N': -54.526406291655, 'O': -74.995458316117,
        'C': -37.788203485235, 'S': -397.663040369707
    },

    'ccsd(t)-f12/cc-pvtz-f12': {
        'H': -0.499946213243, 'N': -54.53000909621, 'O': -75.004127673424,
        'C': -37.789862146471, 'S': -397.675447487865
    },

    'ccsd(t)-f12/cc-pvqz-f12': {
        'H': -0.499994558325, 'N': -54.530515226371, 'O': -75.005600062003,
        'C': -37.789961656228, 'S': -397.676719774973
    },

    'ccsd(t)-f12/cc-pcvdz-f12': {
        'H': -0.499811124128 + SOC['H'], 'N': -54.582137180344 + SOC['N'],
        'O': -75.053045547421 + SOC['O'], 'C': -37.840869118707 + SOC['C']
    },

    'ccsd(t)-f12/cc-pcvtz-f12': {
        'H': -0.499946213243 + SOC['H'], 'N': -54.588545831900 + SOC['N'],
        'O': -75.065995072347 + SOC['O'], 'C': -37.844662139972 + SOC['C']
    },

    'ccsd(t)-f12/cc-pcvqz-f12': {
        'H': -0.499994558325 + SOC['H'], 'N': -54.589137594139 + SOC['N'],
        'O': -75.067412234737 + SOC['O'], 'C': -37.844893820561 + SOC['C']
    },

    'ccsd(t)-f12/cc-pvtz-f12(-pp)': {
        'H': -0.499946213243 + SOC['H'], 'N': -54.53000909621 + SOC['N'],
        'O': -75.004127673424 + SOC['O'], 'C': -37.789862146471 + SOC['C'],
        'S': -397.675447487865 + SOC['S'], 'I': -294.81781766 + SOC['I']
    },

    # ccsd(t)/aug-cc-pvtz(-pp) atomic energies were fit to a set of 8 small molecules:
    # CH4, CH3OH, H2S, H2O, SO2, HI, I2, CH3I
    'ccsd(t)/aug-cc-pvtz(-pp)': {
        'H': -0.499821176024 + SOC['H'], 'O': -74.96738492 + SOC['O'],
        'C': -37.77385697 + SOC['C'], 'S': -397.6461604 + SOC['S'],
        'I': -294.7958443 + SOC['I']
    },

    # note that all atom corrections but S are fitted, the correction for S is calculated
    'ccsd(t)-f12/aug-cc-pvdz': {
        'H': -0.499459066131 + SOC['H'], 'N': -54.524279516472 + SOC['N'],
        'O': -74.992097308083 + SOC['O'], 'C': -37.786694171716 + SOC['C'],
        'S': -397.648733842400 + SOC['S']
    },

    'ccsd(t)-f12/aug-cc-pvtz': {
        'H': -0.499844820798 + SOC['H'], 'N': -54.527419359906 + SOC['N'],
        'O': -75.000001429806 + SOC['O'], 'C': -37.788504810868 + SOC['C'],
        'S': -397.666903000231 + SOC['S']
    },

    'ccsd(t)-f12/aug-cc-pvqz': {
        'H': -0.499949526073 + SOC['H'], 'N': -54.529569719016 + SOC['N'],
        'O': -75.004026586610 + SOC['O'], 'C': -37.789387892348 + SOC['C'],
        'S': -397.671214204994 + SOC['S']
    },

    'b-ccsd(t)-f12/cc-pvdz-f12': {
        'H': -0.499811124128 + SOC['H'], 'N': -54.523269942190 + SOC['N'],
        'O': -74.990725918500 + SOC['O'], 'C': -37.785409916465 + SOC['C'],
        'S': -397.658155086033 + SOC['S']
    },

    'b-ccsd(t)-f12/cc-pvtz-f12': {
        'H': -0.499946213243 + SOC['H'], 'N': -54.528135889213 + SOC['N'],
        'O': -75.001094055506 + SOC['O'], 'C': -37.788233578503 + SOC['C'],
        'S': -397.671745425929 + SOC['S']
    },

    'b-ccsd(t)-f12/cc-pvqz-f12': {
        'H': -0.499994558325 + SOC['H'], 'N': -54.529425753163 + SOC['N'],
        'O': -75.003820485005 + SOC['O'], 'C': -37.789006506290 + SOC['C'],
        'S': -397.674145126931 + SOC['S']
    },

    'b-ccsd(t)-f12/cc-pcvdz-f12': {
        'H': -0.499811124128 + SOC['H'], 'N': -54.578602780288 + SOC['N'],
        'O': -75.048064317367 + SOC['O'], 'C': -37.837592033417 + SOC['C']
    },

    'b-ccsd(t)-f12/cc-pcvtz-f12': {
        'H': -0.499946213243 + SOC['H'], 'N': -54.586402551258 + SOC['N'],
        'O': -75.062767632757 + SOC['O'], 'C': -37.842729156944 + SOC['C']
    },

    'b-ccsd(t)-f12/cc-pcvqz-f12': {
        'H': -0.49999456 + SOC['H'], 'N': -54.587781507581 + SOC['N'],
        'O': -75.065397706471 + SOC['O'], 'C': -37.843634971592 + SOC['C']
    },

    'b-ccsd(t)-f12/aug-cc-pvdz': {
        'H': -0.499459066131 + SOC['H'], 'N': -54.520475581942 + SOC['N'],
        'O': -74.986992215049 + SOC['O'], 'C': -37.783294495799 + SOC['C']
    },

    'b-ccsd(t)-f12/aug-cc-pvtz': {
        'H': -0.499844820798 + SOC['H'], 'N': -54.524927371700 + SOC['N'],
        'O': -74.996328829705 + SOC['O'], 'C': -37.786320700792 + SOC['C']
    },

    'b-ccsd(t)-f12/aug-cc-pvqz': {
        'H': -0.499949526073 + SOC['H'], 'N': -54.528189769291 + SOC['N'],
        'O': -75.001879610563 + SOC['O'], 'C': -37.788165047059 + SOC['C']
    },

    'mp2_rmp2_pvdz': {
        'H': -0.49927840 + SOC['H'], 'N': -54.46141996 + SOC['N'], 'O': -74.89408254 + SOC['O'],
        'C': -37.73792713 + SOC['C']
    },

    'mp2_rmp2_pvtz': {
        'H': -0.49980981 + SOC['H'], 'N': -54.49615972 + SOC['N'], 'O': -74.95506980 + SOC['O'],
        'C': -37.75833104 + SOC['C']
    },

    'mp2_rmp2_pvqz': {
        'H': -0.49994557 + SOC['H'], 'N': -54.50715868 + SOC['N'], 'O': -74.97515364 + SOC['O'],
        'C': -37.76533215 + SOC['C']
    },

    'ccsd-f12/cc-pvdz-f12': {
        'H': -0.499811124128 + SOC['H'], 'N': -54.524325513811 + SOC['N'],
        'O': -74.992326577897 + SOC['O'], 'C': -37.786213495943 + SOC['C']
    },

    'ccsd(t)-f12/cc-pvdz-f12_noscale': {
        'H': -0.499811124128 + SOC['H'], 'N': -54.526026290887 + SOC['N'],
        'O': -74.994751897699 + SOC['O'], 'C': -37.787881871511 + SOC['C']
    },

    'g03_pbepbe_6-311++g_d_p': {
        'H': -0.499812273282 + SOC['H'], 'N': -54.5289567564 + SOC['N'],
        'O': -75.0033596764 + SOC['O'], 'C': -37.7937388736 + SOC['C']
    },

    'fci/cc-pvdz': {
        'C': -37.789527 + SOC['C']
    },

    'fci/cc-pvtz': {
        'C': -37.781266669684 + SOC['C']
    },

    'fci/cc-pvqz': {
        'C': -37.787052110598 + SOC['C']
    },

    # 'bmk/cbsb7' and 'bmk/6-311g(2d,d,p)' have the same corrections
    'bmk/cbsb7': {
        'H': -0.498618853119 + SOC['H'], 'N': -54.5697851544 + SOC['N'],
        'O': -75.0515210278 + SOC['O'], 'C': -37.8287310027 + SOC['C'],
        'P': -341.167615941 + SOC['P'], 'S': -398.001619915 + SOC['S']
    },
    'bmk/6-311g(2d,d,p)': {
        'H': -0.498618853119 + SOC['H'], 'N': -54.5697851544 + SOC['N'],
        'O': -75.0515210278 + SOC['O'], 'C': -37.8287310027 + SOC['C'],
        'P': -341.167615941 + SOC['P'], 'S': -398.001619915 + SOC['S']
    },

    # Fitted to small molecules
    'b3lyp/6-31g(d,p)': {
        'H': -0.500426155, 'C': -37.850331697831, 'O': -75.0535872748806, 'S': -398.100820107242
    },

    # Calculated atomic energies
    'b3lyp/6-311+g(3df,2p)': {
        'H': -0.502155915123 + SOC['H'], 'C': -37.8574709934 + SOC['C'],
        'N': -54.6007233609 + SOC['N'], 'O': -75.0909131284 + SOC['O'],
        'P': -341.281730319 + SOC['P'], 'S': -398.134489850 + SOC['S']
    },

    'wb97x-d/aug-cc-pvtz': {
        'H': -0.502803 + SOC['H'], 'N': -54.585652 + SOC['N'], 'O': -75.068286 + SOC['O'],
        'C': -37.842014 + SOC['C']
    },

    # Calculated atomic energies (unfitted)
    'MRCI+Davidson/aug-cc-pV(T+d)Z': {
        'H': -0.49982118 + SOC['H'], 'C': -37.78321274 + SOC['C'], 'N': -54.51729444 + SOC['N'],
        'O': -74.97847534 + SOC['O'], 'S': -397.6571654 + SOC['S']
    },

}

# Petersson-type bond additivity correction parameters
pbac = {

    # 'S-H', 'C-S', 'C=S', 'S-S', 'O-S', 'O=S', 'O=S=O' taken from http://hdl.handle.net/1721.1/98155 (both for
    # 'CCSD(T)-F12/cc-pVDZ-F12' and 'CCSD(T)-F12/cc-pVTZ-F12')
    'ccsd(t)-f12/cc-pvdz-f12': {
        'C-H': -0.46, 'C-C': -0.68, 'C=C': -1.90, 'C#C': -3.13,
        'O-H': -0.51, 'C-O': -0.23, 'C=O': -0.69, 'O-O': -0.02, 'C-N': -0.67,
        'C=N': -1.46, 'C#N': -2.79, 'N-O': 0.74, 'N_O': -0.23, 'N=O': -0.51,
        'N-H': -0.69, 'N-N': -0.47, 'N=N': -1.54, 'N#N': -2.05, 'S-H': 0.87,
        'C-S': 0.42, 'C=S': 0.51, 'S-S': 0.86, 'O-S': 0.23, 'O=S': -0.53,
        'O=S=O': 1.95
    },

    'ccsd(t)-f12/cc-pvtz-f12': {
        'C-H': -0.09, 'C-C': -0.27, 'C=C': -1.03, 'C#C': -1.79,
        'O-H': -0.06, 'C-O': 0.14, 'C=O': -0.19, 'O-O': 0.16, 'C-N': -0.18,
        'C=N': -0.41, 'C#N': -1.41, 'N-O': 0.87, 'N_O': -0.09, 'N=O': -0.23,
        'N-H': -0.01, 'N-N': -0.21, 'N=N': -0.44, 'N#N': -0.76, 'S-H': 0.52,
        'C-S': 0.13, 'C=S': -0.12, 'S-S': 0.30, 'O-S': 0.15, 'O=S': -2.61,
        'O=S=O': 0.27
    },

    'ccsd(t)-f12/cc-pvqz-f12': {
        'C-H': -0.08, 'C-C': -0.26, 'C=C': -1.01, 'C#C': -1.66,
        'O-H': 0.07, 'C-O': 0.25, 'C=O': -0.03, 'O-O': 0.26, 'C-N': -0.20,
        'C=N': -0.30, 'C#N': -1.33, 'N-O': 1.01, 'N_O': -0.03, 'N=O': -0.26,
        'N-H': 0.06, 'N-N': -0.23, 'N=N': -0.37, 'N#N': -0.64
    },

    'cbs-qb3': {
        'C-H': -0.11, 'C-C': -0.30, 'C=C': -0.08, 'C#C': -0.64, 'O-H': 0.02, 'C-O': 0.33, 'C=O': 0.55,
        # Table IX: Petersson GA (1998) J. of Chemical Physics, DOI: 10.1063/1.477794
        'N-H': -0.42, 'C-N': -0.13, 'C#N': -0.89, 'C-F': 0.55, 'C-Cl': 1.29, 'S-H': 0.0, 'C-S': 0.43,
        'O=S': -0.78,
        'N=O': 1.11, 'N-N': -1.87, 'N=N': -1.58, 'N-O': 0.35,
        # Table 2: Ashcraft R (2007) J. Phys. Chem. B; DOI: 10.1021/jp073539t
        'N#N': -2.0, 'O=O': -0.2, 'H-H': 1.1,  # Unknown source
    },

    'cbs-qb3-paraskevas': {
        # NOTE: The Paraskevas corrections are inaccurate for non-oxygenated hydrocarbons,
        # and may do poorly in combination with the Petersson corrections
        'C-C': -0.495, 'C-H': -0.045, 'C=C': -0.825, 'C-O': 0.378, 'C=O': 0.743, 'O-H': -0.423,
        # Table2: Paraskevas, PD (2013). Chemistry-A European J., DOI: 10.1002/chem.201301381
        'C#C': -0.64, 'C#N': -0.89, 'C-S': 0.43, 'O=S': -0.78, 'S-H': 0.0, 'C-N': -0.13, 'C-Cl': 1.29,
        'C-F': 0.55,  # Table IX: Petersson GA (1998) J. of Chemical Physics, DOI: 10.1063/1.477794
        'N-H': -0.42, 'N=O': 1.11, 'N-N': -1.87, 'N=N': -1.58, 'N-O': 0.35,
        # Table 2: Ashcraft R (2007) J. Phys. Chem. B; DOI: 10.1021/jp073539t
        'N#N': -2.0, 'O=O': -0.2, 'H-H': 1.1,  # Unknown source
    },

    # Identical corrections for 'b3lyp/cbsb7', 'b3lyp/6-311g(2d,d,p)', 'b3lyp/6-311+g(3df,2p)', 'b3lyp/6-31g(d,p)'
    'b3lyp/cbsb7': {
        'C-H': 0.25, 'C-C': -1.89, 'C=C': -0.40, 'C#C': -1.50,
        'O-H': -1.09, 'C-O': -1.18, 'C=O': -0.01, 'N-H': 1.36, 'C-N': -0.44,
        'C#N': 0.22, 'C-S': -2.35, 'O=S': -5.19, 'S-H': -0.52,
    },
    'b3lyp/6-311g(2d,d,p)': {
        'C-H': 0.25, 'C-C': -1.89, 'C=C': -0.40, 'C#C': -1.50,
        'O-H': -1.09, 'C-O': -1.18, 'C=O': -0.01, 'N-H': 1.36, 'C-N': -0.44,
        'C#N': 0.22, 'C-S': -2.35, 'O=S': -5.19, 'S-H': -0.52,
    },
    'b3lyp/6-311+g(3df,2p)': {
        'C-H': 0.25, 'C-C': -1.89, 'C=C': -0.40, 'C#C': -1.50,
        'O-H': -1.09, 'C-O': -1.18, 'C=O': -0.01, 'N-H': 1.36, 'C-N': -0.44,
        'C#N': 0.22, 'C-S': -2.35, 'O=S': -5.19, 'S-H': -0.52,
    },
    'b3lyp/6-31g(d,p)': {
        'C-H': 0.25, 'C-C': -1.89, 'C=C': -0.40, 'C#C': -1.50,
        'O-H': -1.09, 'C-O': -1.18, 'C=O': -0.01, 'N-H': 1.36, 'C-N': -0.44,
        'C#N': 0.22, 'C-S': -2.35, 'O=S': -5.19, 'S-H': -0.52,
    },

    # 'g4' : {
    # 'Br-Br': -4.6523,'Cl-H': -0.8715, 'Cl-O': -0.8703, 'Br-C': -1.884, 'Br-O': -2.1862,
    # 'Br-H': -1.948, 'C-O': 0.5712, 'C-Cl': -0.4162, 'Cl-Cl': -2.0127, 'C-H': 0.7675, 'C-F': 0.0634, 
    # 'C-C': 0.5222, 'C#C': 0.5809, 'O-O': -0.0068, 'C=O': 0.4695, 'H-O': 0.4392, 'H-H': 0.7966, 'C=C': 0.5705,
    # 'F-O': -0.4662, 'F-H': 0.2108, 'F-F': -1.2368, 'O=O': -0.0997
    # },

    # 'g4': {'Br-Br': 0.5635,'Br-C': 0.5283,'Br-Cl': -0.5699,'Br-O': 0.5636,
    # 'C#C': -0.7806,'C#N': -0.4384,'C-C': 0.0824,'C-Cl': -0.6351,'C-F': -0.6744,
    # 'C-H': -0.3155,'C-N': -0.0066,'C-O': 0.1099,'C=C': -0.1923,'C=N': -0.5257,
    # 'C=O': -0.6293,'Cl-Cl': -2.1081,'Cl-F': -1.382,'Cl-H': -1.8078,'Cl-N': -1.6896,
    # 'Cl-O': -1.3419,'F-F': -2.2409,'F-H': -1.1797,'F-N': -0.8146,'F-O': -1.2883,
    # 'H-H': -0.9805,'H-N': -0.6301,'H-O': -0.9369,'N#N': -1.0998,'N-N': 0.0366,
    # 'N-O': -0.0954,'N=N': -0.3661,'N=O': -0.8795,'O-O': -0.7334,'O=O': -1.3686
    # },

    # 0K
    # 'g4' : {'Br-C': -1.1045, 'Br-O': -1.5006, 'C#C': -0.8327,
    #  'C#N': -0.4995, 'C-C': -0.0022, 'C-Cl': -0.6676, 'C-F': -0.6793,
    #  'C-H': -0.3388, 'C-N': 0.0172, 'C-O': 0.16, 'C=C': -0.2126,
    #  'C=N': -0.4442, 'C=O': -0.6359, 'Cl-N': -1.7557, 'Cl-O': -1.3752,
    #  'F-N': -0.8644, 'F-O': -1.3582, 'H-N': -0.6152,
    #  'H-O': -0.8856, 'N#N': -0.9954, 'N-N': 0.0247, 'N-O': -0.1786,
    #  'N=N': -0.2886, 'N=O': -0.8186, 'O-O': -0.8007, 'O=O': -1.3721
    # },

    'g4' : {'Br-C': -2.5224,'Br-O': -2.9547, 'C#C': -0.7397, 'C#N': -0.4711,
        'C-C': 0.1401, 'C-Cl': -0.6857, 'C-F': -0.6821, 'C-H': -0.3766,
        'C-N': 0.0997,'C-O': 0.0638, 'C=C': -0.3035, 'C=N': -0.4056,
        'C=O': -0.5765, 'Cl-N': -1.8766, 'Cl-O': -1.3823, 'F-N': -0.9748,
        'F-O': -1.351, 'H-N': -0.6294, 'H-O': -0.8635, 'N-N': 0.0314,
        'N-O': -0.1821, 'N=N': -0.2762, 'N=O': -0.6745,'O-O': -0.7278,
        'O=O': -4.3844
    },
}

# Melius-type bond additivity correction parameters
mbac = {}
