"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 892 $"


import re

import numpy

import logfileparser
import utils


class GAMESS(logfileparser.Logfile):
    """A GAMESS log file."""
    SCFRMS, SCFMAX, SCFENERGY = range(3) # Used to index self.scftargets[]
    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(GAMESS, self).__init__(logname="GAMESS", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return "GAMESS log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'GAMESS("%s")' % (self.filename)

    def normalisesym(self, label):
        """Normalise the symmetries used by GAMESS.

        To normalise, two rules need to be applied:
        (1) Occurences of U/G in the 2/3 position of the label
            must be lower-cased
        (2) Two single quotation marks must be replaced by a double

        >>> t = GAMESS("dummyfile").normalisesym
        >>> labels = ['A', 'A1', 'A1G', "A'", "A''", "AG"]
        >>> answers = map(t, labels)
        >>> print answers
        ['A', 'A1', 'A1g', "A'", 'A"', 'Ag']
        """
        if label[1:] == "''":
            end = '"'
        else:
            end = label[1:].replace("U", "u").replace("G", "g")
        return label[0] + end

    def before_parsing(self):

        self.firststdorient = True # Used to decide whether to wipe the atomcoords clean
        self.geooptfinished = False # Used to avoid extracting the final geometry twice
        self.cihamtyp = "none" # Type of CI Hamiltonian: saps or dets.
        self.scftype = "none" # Type of SCF calculation: BLYP, RHF, ROHF, etc.
    
    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        if line [1:12] == "INPUT CARD>":
            return

        # We are looking for this line:
        #           PARAMETERS CONTROLLING GEOMETRY SEARCH ARE
        #           ...
        #           OPTTOL = 1.000E-04          RMIN   = 1.500E-03
        if line[10:18] == "OPTTOL =":
            if not hasattr(self, "geotargets"):
                opttol = float(line.split()[2])
                self.geotargets = numpy.array([opttol, 3. / opttol], "d")
                        
        if line.find("FINAL") == 1:
            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
        # Has to deal with such lines as:
        #  FINAL R-B3LYP ENERGY IS     -382.0507446475 AFTER  10 ITERATIONS
        #  FINAL ENERGY IS     -379.7594673378 AFTER   9 ITERATIONS
        # ...so take the number after the "IS"
            temp = line.split()
            self.scfenergies.append(utils.convertor(float(temp[temp.index("IS") + 1]), "hartree", "eV"))

        # Total energies after Moller-Plesset corrections
        if (line.find("RESULTS OF MOLLER-PLESSET") >= 0 or
            line[6:37] == "SCHWARZ INEQUALITY TEST SKIPPED"):
            # Output looks something like this:
            # RESULTS OF MOLLER-PLESSET 2ND ORDER CORRECTION ARE
            #         E(0)=      -285.7568061536
            #         E(1)=         0.0
            #         E(2)=        -0.9679419329
            #       E(MP2)=      -286.7247480864
            # where E(MP2) = E(0) + E(2)
            #
            # with GAMESS-US 12 Jan 2009 (R3) the preceding text is different:
            ##      DIRECT 4-INDEX TRANSFORMATION 
            ##      SCHWARZ INEQUALITY TEST SKIPPED          0 INTEGRAL BLOCKS
            ##                     E(SCF)=       -76.0088477471
            ##                       E(2)=        -0.1403745370
            ##                     E(MP2)=       -76.1492222841            
            if not hasattr(self, "mpenergies"):
                self.mpenergies = []
            # Each iteration has a new print-out
            self.mpenergies.append([])
            # GAMESS-US presently supports only second order corrections (MP2)
            # PC GAMESS also has higher levels (3rd and 4th), with different output
            # Only the highest level MP4 energy is gathered (SDQ or SDTQ)            
            while re.search("DONE WITH MP(\d) ENERGY", line) is None:
                line = inputfile.next()
                if len(line.split()) > 0:
                    # Only up to MP2 correction
                    if line.split()[0] == "E(MP2)=":
                        mp2energy = float(line.split()[1])
                        self.mpenergies[-1].append(utils.convertor(mp2energy, "hartree", "eV"))
                    # MP2 before higher order calculations
                    if line.split()[0] == "E(MP2)":
                        mp2energy = float(line.split()[2])
                        self.mpenergies[-1].append(utils.convertor(mp2energy, "hartree", "eV"))
                    if line.split()[0] == "E(MP3)":
                        mp3energy = float(line.split()[2])
                        self.mpenergies[-1].append(utils.convertor(mp3energy, "hartree", "eV"))
                    if line.split()[0] in ["E(MP4-SDQ)", "E(MP4-SDTQ)"]:
                        mp4energy = float(line.split()[2])
                        self.mpenergies[-1].append(utils.convertor(mp4energy, "hartree", "eV"))

        # Total energies after Coupled Cluster calculations
        # Only the highest Coupled Cluster level result is gathered
        if line[12:23] == "CCD ENERGY:":
            if not hasattr(self, "ccenergies"):
                self.ccenergies = []
            ccenergy = float(line.split()[2])
            self.ccenergies.append(utils.convertor(ccenergy, "hartree", "eV"))
        if line.find("CCSD") >= 0 and line.split()[0:2] == ["CCSD", "ENERGY:"]:
            if not hasattr(self, "ccenergies"):
                self.ccenergies = []
            ccenergy = float(line.split()[2])
            line = inputfile.next()
            if line[8:23] == "CCSD[T] ENERGY:":
                ccenergy = float(line.split()[2])
                line = inputfile.next()
                if line[8:23] == "CCSD(T) ENERGY:":
                    ccenergy = float(line.split()[2])
            self.ccenergies.append(utils.convertor(ccenergy, "hartree", "eV"))
        # Also collect MP2 energies, which are always calculated before CC
        if line [8:23] == "MBPT(2) ENERGY:":
            if not hasattr(self, "mpenergies"):
                self.mpenergies = []
            self.mpenergies.append([])
            mp2energy = float(line.split()[2])
            self.mpenergies[-1].append(utils.convertor(mp2energy, "hartree", "eV"))

        # Extract charge and multiplicity
        if line[1:19] == "CHARGE OF MOLECULE":
            self.charge = int(line.split()[-1])
            self.mult = int(inputfile.next().split()[-1])

        # etenergies (used only for CIS runs now)
        if "EXCITATION ENERGIES" in line and line.find("DONE WITH") < 0:
            if not hasattr(self, "etenergies"):
                self.etenergies = []
            header = inputfile.next().rstrip()
            get_etosc = False
            if header.endswith("OSC. STR."):
                # water_cis_dets.out does not have the oscillator strength
                # in this table...it is extracted from a different section below
                get_etosc = True
                self.etoscs = []
            dashes = inputfile.next()
            line = inputfile.next()
            broken = line.split()
            while len(broken) > 0:
                # Take hartree value with more numbers, and convert.
                # Note that the values listed after this are also less exact!
                etenergy = float(broken[1])
                self.etenergies.append(utils.convertor(etenergy, "hartree", "cm-1"))
                if get_etosc:
                    etosc = float(broken[-1])
                    self.etoscs.append(etosc)
                broken = inputfile.next().split()

        # Detect the CI hamiltonian type, if applicable.
        # Should always be detected if CIS is done.
        if line[8:64] == "RESULTS FROM SPIN-ADAPTED ANTISYMMETRIZED PRODUCT (SAPS)":
            self.cihamtyp = "saps"
        if line[8:64] == "RESULTS FROM DETERMINANT BASED ATOMIC ORBITAL CI-SINGLES":
            self.cihamtyp = "dets"

        # etsecs (used only for CIS runs for now)
        if line[1:14] == "EXCITED STATE":
            if not hasattr(self, 'etsecs'):
                self.etsecs = []
            if not hasattr(self, 'etsyms'):
                self.etsyms = []
            statenumber = int(line.split()[2])
            spin = int(float(line.split()[7]))
            if spin == 0:
                sym = "Singlet"
            if spin == 1:
                sym = "Triplet"
            sym += '-' + line.split()[-1]
            self.etsyms.append(sym)
            # skip 5 lines
            for i in range(5):
                line = inputfile.next()
            line = inputfile.next()
            CIScontribs = []
            while line.strip()[0] != "-":
                MOtype = 0
                # alpha/beta are specified for hamtyp=dets
                if self.cihamtyp == "dets":
                    if line.split()[0] == "BETA":
                        MOtype = 1
                fromMO = int(line.split()[-3])-1
                toMO = int(line.split()[-2])-1
                coeff = float(line.split()[-1])
                # With the SAPS hamiltonian, the coefficients are multiplied
                #   by sqrt(2) so that they normalize to 1.
                # With DETS, both alpha and beta excitations are printed.
                # if self.cihamtyp == "saps":
                #    coeff /= numpy.sqrt(2.0)
                CIScontribs.append([(fromMO,MOtype),(toMO,MOtype),coeff])
                line = inputfile.next()
            self.etsecs.append(CIScontribs)

        # etoscs (used only for CIS runs now)
        if line[1:50] == "TRANSITION FROM THE GROUND STATE TO EXCITED STATE":
            if not hasattr(self, "etoscs"):
                self.etoscs = []
            statenumber = int(line.split()[-1])
            # skip 7 lines
            for i in range(8):
                line = inputfile.next()
            strength = float(line.split()[3])
            self.etoscs.append(strength)

        # TD-DFT for GAMESS-US
        if line[14:29] == "LET EXCITATIONS": # TRIPLET and SINGLET
            self.etenergies = []
            self.etoscs = []
            self.etsecs = []
            etsyms = []
            minus = inputfile.next()
            blank = inputfile.next()
            line = inputfile.next()
            # Loop starts on the STATE line
            while line.find("STATE") >= 0:
                broken = line.split()
                self.etenergies.append(utils.convertor(float(broken[-2]), "eV", "cm-1"))
                broken = inputfile.next().split()
                self.etoscs.append(float(broken[-1]))
                sym = inputfile.next() # Not always present
                if sym.find("SYMMETRY")>=0:
                    etsyms.append(sym.split()[-1])
                    header = inputfile.next()
                minus = inputfile.next()
                CIScontribs = []
                line = inputfile.next()
                while line.strip():
                    broken = line.split()
                    fromMO, toMO = [int(broken[x]) - 1 for x in [2, 4]]
                    CIScontribs.append([(fromMO, 0), (toMO, 0), float(broken[1])])
                    line = inputfile.next()
                self.etsecs.append(CIScontribs)
                line = inputfile.next()
            if etsyms: # Not always present
                self.etsyms = etsyms
         
        # Maximum and RMS gradients.
        if "MAXIMUM GRADIENT" in line or "RMS GRADIENT" in line:

            if not hasattr(self, "geovalues"):
                self.geovalues = []

            parts = line.split()

            # Newer versions (around 2006) have both maximum and RMS on one line:
            #       MAXIMUM GRADIENT =  0.0531540    RMS GRADIENT = 0.0189223
            if len(parts) == 8:
                maximum = float(parts[3])
                rms = float(parts[7])
            
            # In older versions of GAMESS, this spanned two lines, like this:
            #       MAXIMUM GRADIENT =    0.057578167
            #           RMS GRADIENT =    0.027589766
            if len(parts) == 4:
                maximum = float(parts[3])
                line = inputfile.next()
                parts = line.split()
                rms = float(parts[3])


            # FMO also prints two final one- and two-body gradients (see exam37):
            #   (1) MAXIMUM GRADIENT =  0.0531540    RMS GRADIENT = 0.0189223
            if len(parts) == 9:
                maximum = float(parts[4])
                rms = float(parts[8])

            self.geovalues.append([maximum, rms])

        if line[11:50] == "ATOMIC                      COORDINATES":
            # This is the input orientation, which is the only data available for
            # SP calcs, but which should be overwritten by the standard orientation
            # values, which is the only information available for all geoopt cycles.
            if not hasattr(self, "atomcoords"):
                self.atomcoords = []
                self.atomnos = []
            line = inputfile.next()
            atomcoords = []
            atomnos = []
            line = inputfile.next()
            while line.strip():
                temp = line.strip().split()
                atomcoords.append([utils.convertor(float(x), "bohr", "Angstrom") for x in temp[2:5]])
                atomnos.append(int(round(float(temp[1])))) # Don't use the atom name as this is arbitary
                line = inputfile.next()
            self.atomnos = numpy.array(atomnos, "i")
            self.atomcoords.append(atomcoords)

        if line[12:40] == "EQUILIBRIUM GEOMETRY LOCATED":
            # Prevent extraction of the final geometry twice
            self.geooptfinished = True
        
        if line[1:29] == "COORDINATES OF ALL ATOMS ARE" and not self.geooptfinished:
            # This is the standard orientation, which is the only coordinate
            # information available for all geometry optimisation cycles.
            # The input orientation will be overwritten if this is a geometry optimisation
            # We assume that a previous Input Orientation has been found and
            # used to extract the atomnos
            if self.firststdorient:
                self.firststdorient = False
                # Wipes out the single input coordinate at the start of the file
                self.atomcoords = []
                
            line = inputfile.next()
            hyphens = inputfile.next()

            atomcoords = []
            line = inputfile.next()                

            for i in range(self.natom):
                temp = line.strip().split()
                atomcoords.append(map(float, temp[2:5]))
                line = inputfile.next()
            self.atomcoords.append(atomcoords)
        
        # Section with SCF information.
        #
        # The space at the start of the search string is to differentiate from MCSCF.
        # Everything before the search string is stored as the type of SCF.
        # SCF types may include: BLYP, RHF, ROHF, UHF, etc.
        #
        # For example, in exam17 the section looks like this (note that this is GVB):
        #          ------------------------
        #          ROHF-GVB SCF CALCULATION
        #          ------------------------
        # GVB STEP WILL USE    119875 WORDS OF MEMORY.
        #
        #     MAXIT=  30   NPUNCH= 2   SQCDF TOL=1.0000E-05
        #     NUCLEAR ENERGY=        6.1597411978
        #     EXTRAP=T   DAMP=F   SHIFT=F   RSTRCT=F   DIIS=F  SOSCF=F
        #
        # ITER EX     TOTAL ENERGY       E CHANGE        SQCDF       DIIS ERROR
        #   0  0      -38.298939963   -38.298939963   0.131784454   0.000000000
        #   1  1      -38.332044339    -0.033104376   0.026019716   0.000000000
        # ... and will be terminated by a blank line.
        if line.rstrip()[-16:] == " SCF CALCULATION":

            # Remember the type of SCF.
            self.scftype = line.strip()[:-16]

            dashes = inputfile.next()

            while line [:5] != " ITER":

                # GVB uses SQCDF for checking convergence (for example in exam17).
                if "GVB" in self.scftype and "SQCDF TOL=" in line:
                    scftarget = float(line.split("=")[-1])

                # Normally however the density is used as the convergence criterium.
                # Deal with various versions:
                #   (GAMESS VERSION = 12 DEC 2003)
                #     DENSITY MATRIX CONV=  2.00E-05  DFT GRID SWITCH THRESHOLD=  3.00E-04
                #   (GAMESS VERSION = 22 FEB 2006)
                #     DENSITY MATRIX CONV=  1.00E-05
                #   (PC GAMESS version 6.2, Not DFT?)
                #     DENSITY CONV=  1.00E-05
                elif "DENSITY CONV" in line or "DENSITY MATRIX CONV" in line:
                    scftarget = float(line.split()[-1])

                line = inputfile.next()

            if not hasattr(self, "scftargets"):
                self.scftargets = []

            self.scftargets.append([scftarget])

            if not hasattr(self,"scfvalues"):
                self.scfvalues = []

            line = inputfile.next()

            # Normally the iteration print in 6 columns.
            # For ROHF, however, it is 5 columns, thus this extra parameter.
            if "ROHF" in self.scftype:
                valcol = 4
            else:
                valcol = 5

            # SCF iterations are terminated by a blank line.
            # The first four characters usually contains the step number.
            # However, lines can also contain messages, including:
            #   * * *   INITIATING DIIS PROCEDURE   * * *
            #   CONVERGED TO SWOFF, SO DFT CALCULATION IS NOW SWITCHED ON
            #   DFT CODE IS SWITCHING BACK TO THE FINER GRID
            values = []
            while line.strip():
                try:
                    temp = int(line[0:4])
                except ValueError:
                    pass
                else:
                    values.append([float(line.split()[valcol])])
                line = inputfile.next()
            self.scfvalues.append(values)

        if line.find("NORMAL COORDINATE ANALYSIS IN THE HARMONIC APPROXIMATION") >= 0:
        # GAMESS has...
        # MODES 1 TO 6 ARE TAKEN AS ROTATIONS AND TRANSLATIONS.
        #
        #     FREQUENCIES IN CM**-1, IR INTENSITIES IN DEBYE**2/AMU-ANGSTROM**2,
        #     REDUCED MASSES IN AMU.
        #
        #                          1           2           3           4           5
        #       FREQUENCY:        52.49       41.45       17.61        9.23       10.61  
        #    REDUCED MASS:      3.92418     3.77048     5.43419     6.44636     5.50693
        #    IR INTENSITY:      0.00013     0.00001     0.00004     0.00000     0.00003

        # ...or in the case of a numerical Hessian job...

        # MODES 1 TO 5 ARE TAKEN AS ROTATIONS AND TRANSLATIONS.
        #
        #     FREQUENCIES IN CM**-1, IR INTENSITIES IN DEBYE**2/AMU-ANGSTROM**2,
        #     REDUCED MASSES IN AMU.
        #
        #                          1           2           3           4           5
        #       FREQUENCY:         0.05        0.03        0.03       30.89       30.94  
        #    REDUCED MASS:      8.50125     8.50137     8.50136     1.06709     1.06709

        
        # whereas PC-GAMESS has...
        # MODES 1 TO 6 ARE TAKEN AS ROTATIONS AND TRANSLATIONS.
        #
        #     FREQUENCIES IN CM**-1, IR INTENSITIES IN DEBYE**2/AMU-ANGSTROM**2
        #
        #                          1           2           3           4           5
        #       FREQUENCY:         5.89        1.46        0.01        0.01        0.01  
        #    IR INTENSITY:      0.00000     0.00000     0.00000     0.00000     0.00000
        
        # If Raman is present we have (for PC-GAMESS)...
        # MODES 1 TO 6 ARE TAKEN AS ROTATIONS AND TRANSLATIONS.
        #
        #     FREQUENCIES IN CM**-1, IR INTENSITIES IN DEBYE**2/AMU-ANGSTROM**2
        #     RAMAN INTENSITIES IN ANGSTROM**4/AMU, DEPOLARIZATIONS ARE DIMENSIONLESS
        #
        #                          1           2           3           4           5
        #       FREQUENCY:         5.89        1.46        0.04        0.03        0.01  
        #    IR INTENSITY:      0.00000     0.00000     0.00000     0.00000     0.00000
        # RAMAN INTENSITY:       12.675       1.828       0.000       0.000       0.000
        #  DEPOLARIZATION:        0.750       0.750       0.124       0.009       0.750

        # If PC-GAMESS has not reached the stationary point we have
        # MODES 1 TO 5 ARE TAKEN AS ROTATIONS AND TRANSLATIONS.
        #
        #     FREQUENCIES IN CM**-1, IR INTENSITIES IN DEBYE**2/AMU-ANGSTROM**2
        #
        #     *******************************************************
        #     * THIS IS NOT A STATIONARY POINT ON THE MOLECULAR PES *
        #     *     THE VIBRATIONAL ANALYSIS IS NOT VALID !!!       *
        #     *******************************************************
        #
        #                          1           2           3           4           5
        
        # MODES 2 TO 7 ARE TAKEN AS ROTATIONS AND TRANSLATIONS.

            self.vibfreqs = []
            self.vibirs = []
            self.vibdisps = []

            # Need to get to the modes line
            warning = False
            while line.find("MODES") == -1:
                line = inputfile.next()
                if line.find("THIS IS NOT A STATIONARY POINT")>=0:
                    warning = True
            startrot = int(line.split()[1])
            endrot = int(line.split()[3])
            blank = inputfile.next()

            line = inputfile.next() # FREQUENCIES, etc.
            while line != blank:
                line = inputfile.next()
            if warning: # Get past the second warning
                line = inputfile.next()
                while line!= blank:
                    line = inputfile.next()
                self.logger.warning("This is not a stationary point on the molecular"
                                    "PES. The vibrational analysis is not valid.")
            
            freqNo = inputfile.next()
            while freqNo.find("SAYVETZ") == -1:
                freq = inputfile.next().strip().split()[1:]
            # May include imaginary frequencies
            #       FREQUENCY:       825.18 I    111.53       12.62       10.70        0.89
                newfreq = []
                for i, x in enumerate(freq):
                    if x!="I":
                        newfreq.append(float(x))
                    else:
                        newfreq[-1] = -newfreq[-1]
                self.vibfreqs.extend(newfreq)
                line = inputfile.next()
                if line.find("REDUCED") >= 0: # skip the reduced mass (not always present)
                    line = inputfile.next()
                if line.find("IR INTENSITY") >= 0:
                    # Not present if a numerical Hessian calculation
                    irIntensity = map(float, line.strip().split()[2:])
                    self.vibirs.extend([utils.convertor(x, "Debye^2/amu-Angstrom^2", "km/mol") for x in irIntensity])
                    line = inputfile.next()
                if line.find("RAMAN") >= 0:
                    if not hasattr(self,"vibramans"):
                        self.vibramans = []
                    ramanIntensity = line.strip().split()
                    self.vibramans.extend(map(float, ramanIntensity[2:]))
                    depolar = inputfile.next()
                    line = inputfile.next()
                assert line == blank

                # Extract the Cartesian displacement vectors
                p = [ [], [], [], [], [] ]
                for j in range(len(self.atomnos)):
                    q = [ [], [], [], [], [] ]
                    for k in range(3): # x, y, z
                        line = inputfile.next()[21:]
                        broken = map(float, line.split())
                        for l in range(len(broken)):
                            q[l].append(broken[l])
                    for k in range(len(broken)):
                        p[k].append(q[k])
                self.vibdisps.extend(p[:len(broken)])

                # Skip the Sayvetz stuff at the end
                for j in range(10):
                    line = inputfile.next()
                blank = inputfile.next()
                freqNo = inputfile.next()
            # Exclude rotations and translations
            self.vibfreqs = numpy.array(self.vibfreqs[:startrot-1]+self.vibfreqs[endrot:], "d")
            self.vibirs = numpy.array(self.vibirs[:startrot-1]+self.vibirs[endrot:], "d")
            self.vibdisps = numpy.array(self.vibdisps[:startrot-1]+self.vibdisps[endrot:], "d")
            if hasattr(self, "vibramans"):
                self.vibramans = numpy.array(self.vibramans[:startrot-1]+self.vibramans[endrot:], "d")

        if line[5:21] == "ATOMIC BASIS SET":
            self.gbasis = []
            line = inputfile.next()
            while line.find("SHELL")<0:
                line = inputfile.next()
            blank = inputfile.next()
            atomname = inputfile.next()
            # shellcounter stores the shell no of the last shell
            # in the previous set of primitives
            shellcounter = 1
            while line.find("TOTAL NUMBER")<0:
                blank = inputfile.next()
                line = inputfile.next()
                shellno = int(line.split()[0])
                shellgap = shellno - shellcounter
                gbasis = [] # Stores basis sets on one atom
                shellsize = 0
                while len(line.split())!=1 and line.find("TOTAL NUMBER")<0:
                    shellsize += 1
                    coeff = {}
                    # coefficients and symmetries for a block of rows
                    while line.strip():
                        temp = line.strip().split()
                        sym = temp[1]
                        assert sym in ['S', 'P', 'D', 'F', 'G', 'L']
                        if sym == "L": # L refers to SP
                            if len(temp)==6: # GAMESS US
                                coeff.setdefault("S", []).append( (float(temp[3]), float(temp[4])) )
                                coeff.setdefault("P", []).append( (float(temp[3]), float(temp[5])) )
                            else: # PC GAMESS
                                assert temp[6][-1] == temp[9][-1] == ')'
                                coeff.setdefault("S", []).append( (float(temp[3]), float(temp[6][:-1])) )
                                coeff.setdefault("P", []).append( (float(temp[3]), float(temp[9][:-1])) )
                        else:
                            if len(temp)==5: # GAMESS US
                                coeff.setdefault(sym, []).append( (float(temp[3]), float(temp[4])) )
                            else: # PC GAMESS
                                assert temp[6][-1] == ')'
                                coeff.setdefault(sym, []).append( (float(temp[3]), float(temp[6][:-1])) )
                        line = inputfile.next()
                    # either a blank or a continuation of the block
                    if sym == "L":
                        gbasis.append( ('S', coeff['S']))
                        gbasis.append( ('P', coeff['P']))
                    else:
                        gbasis.append( (sym, coeff[sym]))
                    line = inputfile.next()
                # either the start of the next block or the start of a new atom or
                # the end of the basis function section
                
                numtoadd = 1 + (shellgap / shellsize)
                shellcounter = shellno + shellsize
                for x in range(numtoadd):
                    self.gbasis.append(gbasis)

        if line.find("EIGENVECTORS") == 10 or line.find("MOLECULAR OBRITALS") == 10:
            # The details returned come from the *final* report of evalues and
            #   the last list of symmetries in the log file.
            # Should be followed by lines like this:
            #           ------------
            #           EIGENVECTORS
            #           ------------
            # 
            #                       1          2          3          4          5
            #                   -10.0162   -10.0161   -10.0039   -10.0039   -10.0029
            #                      BU         AG         BU         AG         AG  
            #     1  C  1  S    0.699293   0.699290  -0.027566   0.027799   0.002412
            #     2  C  1  S    0.031569   0.031361   0.004097  -0.004054  -0.000605
            #     3  C  1  X    0.000908   0.000632  -0.004163   0.004132   0.000619
            #     4  C  1  Y   -0.000019   0.000033   0.000668  -0.000651   0.005256
            #     5  C  1  Z    0.000000   0.000000   0.000000   0.000000   0.000000
            #     6  C  2  S   -0.699293   0.699290   0.027566   0.027799   0.002412
            #     7  C  2  S   -0.031569   0.031361  -0.004097  -0.004054  -0.000605
            #     8  C  2  X    0.000908  -0.000632  -0.004163  -0.004132  -0.000619
            #     9  C  2  Y   -0.000019  -0.000033   0.000668   0.000651  -0.005256
            #    10  C  2  Z    0.000000   0.000000   0.000000   0.000000   0.000000
            #    11  C  3  S   -0.018967  -0.019439   0.011799  -0.014884  -0.452328
            #    12  C  3  S   -0.007748  -0.006932   0.000680  -0.000695  -0.024917
            #    13  C  3  X    0.002628   0.002997   0.000018   0.000061  -0.003608
            # and so forth... with blanks lines between blocks of 5 orbitals each.
            # Warning! There are subtle differences between GAMESS-US and PC-GAMES
            #   in the formatting of the first four columns.
            #
            # Watch out for F orbitals...
            # PC GAMESS
            #   19  C   1 YZ   0.000000   0.000000   0.000000   0.000000   0.000000
            #   20  C    XXX   0.000000   0.000000   0.000000   0.000000   0.002249
            #   21  C    YYY   0.000000   0.000000  -0.025555   0.000000   0.000000
            #   22  C    ZZZ   0.000000   0.000000   0.000000   0.002249   0.000000
            #   23  C    XXY   0.000000   0.000000   0.001343   0.000000   0.000000
            # GAMESS US
            #   55  C  1 XYZ   0.000000   0.000000   0.000000   0.000000   0.000000
            #   56  C  1XXXX  -0.000014  -0.000067   0.000000   0.000000   0.000000
            #
            # This is fine for GeoOpt and SP, but may be weird for TD and Freq.

            # This is the stuff that we can read from these blocks.
            self.moenergies = [[]]
            self.mosyms = [[]]
            if not hasattr(self, "nmo"):
                self.nmo = self.nbasis
            self.mocoeffs = [numpy.zeros((self.nmo, self.nbasis), "d")]
            readatombasis = False
            if not hasattr(self, "atombasis"):
                self.atombasis = []
                self.aonames = []
                for i in range(self.natom):
                    self.atombasis.append([])
                self.aonames = []
                readatombasis = True

            dashes = inputfile.next()
            for base in range(0, self.nmo, 5):

                line = inputfile.next()
                # Make sure that this section does not end prematurely - checked by regression test 2CO.ccsd.aug-cc-pVDZ.out.
                if line.strip() != "":
                    break;
                
                numbers = inputfile.next() # Eigenvector numbers.

                # Sometimes there are some blank lines here.
                while not line.strip():
                    line = inputfile.next()

                # Eigenvalues for these orbitals (in hartrees).
                try:
                    self.moenergies[0].extend([utils.convertor(float(x), "hartree", "eV") for x in line.split()])
                except:
                    self.logger.warning('MO section found but could not be parsed!')
                    break;

                # Orbital symmetries.
                line = inputfile.next()
                if line.strip():
                    self.mosyms[0].extend(map(self.normalisesym, line.split()))
                
                # Now we have nbasis lines.
                # Going to use the same method as for normalise_aonames()
                # to extract basis set information.
                p = re.compile("(\d+)\s*([A-Z][A-Z]?)\s*(\d+)\s*([A-Z]+)")
                oldatom ='0'
                for i in range(self.nbasis):
                    line = inputfile.next()

                    # If line is empty, break (ex. for FMO in exam37).
                    if not line.strip(): break

                    # Fill atombasis and aonames only first time around
                    if readatombasis and base == 0:
                        aonames = []
                        start = line[:17].strip()
                        m = p.search(start)
                        if m:
                            g = m.groups()
                            aoname = "%s%s_%s" % (g[1].capitalize(), g[2], g[3])
                            oldatom = g[2]
                            atomno = int(g[2])-1
                            orbno = int(g[0])-1
                        else: # For F orbitals, as shown above
                            g = [x.strip() for x in line.split()]
                            aoname = "%s%s_%s" % (g[1].capitalize(), oldatom, g[2])
                            atomno = int(oldatom)-1
                            orbno = int(g[0])-1
                        self.atombasis[atomno].append(orbno)
                        self.aonames.append(aoname)
                    coeffs = line[15:] # Strip off the crud at the start.
                    j = 0
                    while j*11+4 < len(coeffs):
                        self.mocoeffs[0][base+j, i] = float(coeffs[j * 11:(j + 1) * 11])
                        j += 1

            line = inputfile.next()
            # If it's restricted and no more properties:
            #  ...... END OF RHF/DFT CALCULATION ......
            # If there are more properties (DENSITY MATRIX):
            #               --------------
            #
            # If it's unrestricted we have:
            #
            #  ----- BETA SET ----- 
            #
            #          ------------
            #          EIGENVECTORS
            #          ------------
            #
            #                      1          2          3          4          5
            # ... and so forth.
            line = inputfile.next()
            if line[2:22] == "----- BETA SET -----":
                self.mocoeffs.append(numpy.zeros((self.nmo, self.nbasis), "d"))
                self.moenergies.append([])
                self.mosyms.append([])
                for i in range(4):
                    line = inputfile.next()
                for base in range(0, self.nmo, 5):
                    blank = inputfile.next()
                    line = inputfile.next() # Eigenvector no
                    line = inputfile.next()
                    self.moenergies[1].extend([utils.convertor(float(x), "hartree", "eV") for x in line.split()])
                    line = inputfile.next()
                    self.mosyms[1].extend(map(self.normalisesym, line.split()))
                    for i in range(self.nbasis):
                        line = inputfile.next()
                        temp = line[15:] # Strip off the crud at the start
                        j = 0
                        while j * 11 + 4 < len(temp):
                            self.mocoeffs[1][base+j, i] = float(temp[j * 11:(j + 1) * 11])
                            j += 1
                line = inputfile.next()
            self.moenergies = [numpy.array(x, "d") for x in self.moenergies]

        # Natural orbitals - presently support only CIS.
        # Looks basically the same as eigenvectors, without symmetry labels.
        if line[10:30] == "CIS NATURAL ORBITALS":

            self.nocoeffs = numpy.zeros((self.nmo, self.nbasis), "d")

            dashes = inputfile.next()
            for base in range(0, self.nmo, 5):

                blank = inputfile.next()
                numbers = inputfile.next() # Eigenvector numbers.

                # Eigenvalues for these natural orbitals (not in hartrees!).
                # Sometimes there are some blank lines before it.
                line = inputfile.next()
                while not line.strip():
                    line = inputfile.next()
                eigenvalues = line

                # Orbital symemtry labels are normally here for MO coefficients.
                line = inputfile.next()
                
                # Now we have nbasis lines with the coefficients.
                for i in range(self.nbasis):

                    line = inputfile.next()
                    coeffs = line[15:]
                    j = 0
                    while j*11+4 < len(coeffs):
                        self.nocoeffs[base+j, i] = float(coeffs[j * 11:(j + 1) * 11])
                        j += 1

        # We cannot trust this self.homos until we come to the phrase:
        #   SYMMETRIES FOR INITAL GUESS ORBITALS FOLLOW
        # which either is followed by "ALPHA" or "BOTH" at which point we can say
        # for certain that it is an un/restricted calculations.
        # Note that MCSCF calcs also print this search string, so make sure
        #   that self.homos does not exist yet.
        if line[1:28] == "NUMBER OF OCCUPIED ORBITALS" and not hasattr(self,'homos'):
            homos = [int(line.split()[-1])-1]
            line = inputfile.next()
            homos.append(int(line.split()[-1])-1)
            self.homos = numpy.array(homos, "i")

        
        if line.find("SYMMETRIES FOR INITIAL GUESS ORBITALS FOLLOW") >= 0:
            # Not unrestricted, so lop off the second index.
            # In case the search string above was not used (ex. FMO in exam38),
            #   we can try to use the next line which should also contain the
            #   number of occupied orbitals.
            if line.find("BOTH SET(S)") >= 0:
                nextline = inputfile.next()
                if "ORBITALS ARE OCCUPIED" in nextline:
                    homos = int(nextline.split()[0])-1
                    if hasattr(self,"homos"):
                        try:
                            assert self.homos[0] == homos
                        except AssertionError:
                            self.logger.warning("Number of occupied orbitals not consistent. This is normal for ECP and FMO jobs.")
                    else:
                        self.homos = [homos]
                self.homos = numpy.resize(self.homos, [1])

        # Set the total number of atoms, only once.
        # Normally GAMESS print TOTAL NUMBER OF ATOMS, however in some cases
        #   this is slightly different (ex. lower case for FMO in exam37).
        if not hasattr(self,"natom") and "NUMBER OF ATOMS" in line.upper():
            self.natom = int(line.split()[-1])
            
        if line.find("NUMBER OF CARTESIAN GAUSSIAN BASIS") == 1 or line.find("TOTAL NUMBER OF BASIS FUNCTIONS") == 1:
            # The first is from Julien's Example and the second is from Alexander's
            # I think it happens if you use a polar basis function instead of a cartesian one
            self.nbasis = int(line.strip().split()[-1])
                
        elif line.find("SPHERICAL HARMONICS KEPT IN THE VARIATION SPACE") >= 0:
            # Note that this line is present if ISPHER=1, e.g. for C_bigbasis
            self.nmo = int(line.strip().split()[-1])
            
        elif line.find("TOTAL NUMBER OF MOS IN VARIATION SPACE") == 1:
            # Note that this line is not always present, so by default
            # NBsUse is set equal to NBasis (see below).
            self.nmo = int(line.split()[-1])

        elif line.find("OVERLAP MATRIX") == 0 or line.find("OVERLAP MATRIX") == 1:
            # The first is for PC-GAMESS, the second for GAMESS
            # Read 1-electron overlap matrix
            if not hasattr(self, "aooverlaps"):
                self.aooverlaps = numpy.zeros((self.nbasis, self.nbasis), "d")
            else:
                self.logger.info("Reading additional aooverlaps...")
            base = 0
            while base < self.nbasis:
                blank = inputfile.next()
                line = inputfile.next() # Basis fn number
                blank = inputfile.next()
                for i in range(self.nbasis - base): # Fewer lines each time
                    line = inputfile.next()
                    temp = line.split()
                    for j in range(4, len(temp)):
                        self.aooverlaps[base+j-4, i+base] = float(temp[j])
                        self.aooverlaps[i+base, base+j-4] = float(temp[j])
                base += 5

        # ECP Pseudopotential information
        if "ECP POTENTIALS" in line:
            if not hasattr(self, "coreelectrons"):
                self.coreelectrons = [0]*self.natom
            dashes = inputfile.next()
            blank = inputfile.next()
            header = inputfile.next()
            while header.split()[0] == "PARAMETERS":
                name = header[17:25]
                atomnum = int(header[34:40])
                # The pseudopotnetial is given explicitely
                if header[40:50] == "WITH ZCORE":
                  zcore = int(header[50:55])
                  lmax = int(header[63:67])
                  self.coreelectrons[atomnum-1] = zcore
                # The pseudopotnetial is copied from another atom
                if header[40:55] == "ARE THE SAME AS":
                  atomcopy = int(header[60:])
                  self.coreelectrons[atomnum-1] = self.coreelectrons[atomcopy-1]
                line = inputfile.next()
                while line.split() <> []:
                    line = inputfile.next()
                header = inputfile.next()

        # This was used before refactoring the parser, geotargets was set here after parsing.
        #if not hasattr(self, "geotargets"):
        #    opttol = 1e-4
        #    self.geotargets = numpy.array([opttol, 3. / opttol], "d")
        #if hasattr(self,"geovalues"): self.geovalues = numpy.array(self.geovalues, "d")

        
if __name__ == "__main__":
    import doctest, gamessparser
    doctest.testmod(gamessparser, verbose=False)
