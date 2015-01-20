# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Parser for Q-Chem output files"""

from __future__ import print_function

import re
import numpy

from . import logfileparser
from . import utils


class QChem(logfileparser.Logfile):
    """A Q-Chem 4 log file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(QChem, self).__init__(logname="QChem", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return "QChem log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'QChem("%s")' % (self.filename)

    def normalisesym(self, label):
        """Q-Chem does not require normalizing symmetry labels."""

    def before_parsing(self):

        # Keep track of whether or not we're performing an
        # (un)restricted calculation.
        self.unrestricted = False

        # Compile the dashes-and-or-spaces-only regex.
        self.dashes_and_spaces = re.compile('^[\s-]+$')

    def after_parsing(self):

        # If parsing a fragment job, each of the geometries appended to
        # `atomcoords` may be of different lengths, which will prevent
        # conversion from a list to NumPy array.
        # Take the length of the first geometry as correct, and remove
        # all others with different lengths.
        if len(self.atomcoords) > 1:
            correctlen = len(self.atomcoords[0])
            self.atomcoords[:] = [coords for coords in self.atomcoords
                                  if len(coords) == correctlen]
        # At the moment, there is no similar correction for other properties!

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        # Charge and multiplicity.
        # Only present in the input file.

        if '$molecule' in line:
            line = next(inputfile)
            charge, mult = map(int, line.split())
            self.set_attribute('charge', charge)
            self.set_attribute('mult', mult)

        # Extract the atomic numbers and coordinates of the atoms.

        if 'Standard Nuclear Orientation (Angstroms)' in line:
            if not hasattr(self, 'atomcoords'):
                self.atomcoords = []
            self.skip_lines(inputfile, ['cols', 'dashes'])
            atomelements = []
            atomcoords = []
            line = next(inputfile)
            while list(set(line.strip())) != ['-']:
                entry = line.split()
                atomelements.append(entry[1])
                atomcoords.append(list(map(float, entry[2:])))
                line = next(inputfile)

            self.atomcoords.append(atomcoords)

            if not hasattr(self, 'atomnos'):
                self.atomnos = []
                for atomelement in atomelements:
                    if atomelement == 'GH':
                        self.atomnos.append(0)
                    else:
                        self.atomnos.append(utils.PeriodicTable().number[atomelement])
                self.natom = len(self.atomnos)

        # Number of electrons.
        # Useful for determining the number of occupied/virtual orbitals.
        if 'Nuclear Repulsion Energy' in line:
            if not hasattr(self, 'nalpha'):
                line = next(inputfile)
                nelec_re_string = 'There are(\s+[0-9]+) alpha and(\s+[0-9]+) beta electrons'
                match = re.findall(nelec_re_string, line.strip())
                self.nalpha = int(match[0][0].strip())
                self.nbeta = int(match[0][1].strip())

        # Number of basis functions.
        # Because Q-Chem's integral recursion scheme is defined using
        # Cartesian basis functions, there is often a distinction between the
        # two in the output. We only parse for *pure* functions.
        # Examples:
        #  Only one type:
        #   There are 30 shells and 60 basis functions
        #  Both Cartesian and pure:
        #   ...

        if 'basis functions' in line:
            self.set_attribute('nbasis', int(line.split()[-3]))

        # Check for whether or not we're peforming an
        # (un)restricted calculation.
        if 'calculation will be' in line:
            if ' restricted' in line:
                self.unrestricted = False
            if 'unrestricted' in line:
                self.unrestricted = True

        # Section with SCF iterations goes like this:
        #
        # SCF converges when DIIS error is below 1.0E-05
        # ---------------------------------------
        #  Cycle       Energy         DIIS Error
        # ---------------------------------------
        #    1    -381.9238072190      1.39E-01
        #    2    -382.2937212775      3.10E-03
        #    3    -382.2939780242      3.37E-03
        # ...
        #
        if 'SCF converges when ' in line:
            if not hasattr(self, 'scftargets'):
                self.scftargets = []
            target = float(line.split()[-1])
            self.scftargets.append([target])

            # We should have the header between dashes now,
            # but sometimes there are lines before the first dashes.
            while not 'Cycle       Energy' in line:
                line = next(inputfile)
            self.skip_lines(inputfile, ['d'])

            values = []
            iter_counter = 1
            line = next(inputfile)
            while 'energy in the final basis set' not in line:

                # Some trickery to avoid a lot of printing that can occur
                # between each SCF iteration.
                entry = line.split()
                if len(entry) > 0:
                    if entry[0] == str(iter_counter):
                        # Q-Chem only outputs one error metric.
                        error = float(entry[2])
                        values.append([error])
                        iter_counter += 1
                line = next(inputfile)

                # This is printed in regression Qchem4.2/dvb_sp_unconverged.out
                # so use it to bail out when convergence fails.
                if "SCF failed to converge" in line or "Convergence failure" in line:
                    break

            if not hasattr(self, 'scfvalues'):
                self.scfvalues = []
            self.scfvalues.append(numpy.array(values))

        if 'Total energy in the final basis set' in line:
            if not hasattr(self, 'scfenergies'):
                self.scfenergies = []
            scfenergy = float(line.split()[-1])
            self.scfenergies.append(utils.convertor(scfenergy, 'hartree', 'eV'))

        # Geometry optimization.

        if 'Maximum     Tolerance    Cnvgd?' in line:
            line_g = list(map(float, next(inputfile).split()[1:3]))
            line_d = list(map(float, next(inputfile).split()[1:3]))
            line_e = next(inputfile).split()[2:4]

            if not hasattr(self, 'geotargets'):
                self.geotargets = [line_g[1], line_d[1], self.float(line_e[1])]
            if not hasattr(self, 'geovalues'):
                self.geovalues = []
            try:
                ediff = abs(self.float(line_e[0]))
            except ValueError:
                ediff = numpy.nan
            geovalues = [line_g[0], line_d[0], ediff]
            self.geovalues.append(geovalues)

        if '**  OPTIMIZATION CONVERGED  **' in line:
            if not hasattr(self, 'optdone'):
                self.optdone = []
            self.optdone.append(len(self.atomcoords))

        if '**  MAXIMUM OPTIMIZATION CYCLES REACHED  **' in line:
            if not hasattr(self, 'optdone'):
                self.optdone = []

        # Moller-Plesset corrections.

        # There are multiple modules in Q-Chem for calculating MPn energies:
        # cdman, ccman, and ccman2, all with different output.
        #
        # MP2, RI-MP2, and local MP2 all default to cdman, which has a simple
        # block of output after the regular SCF iterations.
        #
        # MP3 is handled by ccman2.
        #
        # MP4 and variants are handled by ccman.

        # This is the MP2/cdman case.
        if 'MP2         total energy' in line:
            if not hasattr(self, 'mpenergies'):
                self.mpenergies = []
            mp2energy = float(line.split()[4])
            mp2energy = utils.convertor(mp2energy, 'hartree', 'eV')
            self.mpenergies.append([mp2energy])

        # This is the MP3/ccman2 case.
        if line[1:11] == 'MP2 energy' and line[12:19] != 'read as':
            if not hasattr(self, 'mpenergies'):
                self.mpenergies = []
            mpenergies = []
            mp2energy = float(line.split()[3])
            mpenergies.append(mp2energy)
            line = next(inputfile)
            line = next(inputfile)
            # Just a safe check.
            if 'MP3 energy' in line:
                mp3energy = float(line.split()[3])
                mpenergies.append(mp3energy)
            mpenergies = [utils.convertor(mpe, 'hartree', 'eV')
                          for mpe in mpenergies]
            self.mpenergies.append(mpenergies)

        # This is the MP4/ccman case.
        if 'EHF' in line:
            if not hasattr(self, 'mpenergies'):
                self.mpenergies = []
            mpenergies = []

            while list(set(line.strip())) != ['-']:

                if 'EMP2' in line:
                    mp2energy = float(line.split()[2])
                    mpenergies.append(mp2energy)
                if 'EMP3' in line:
                    mp3energy = float(line.split()[2])
                    mpenergies.append(mp3energy)
                if 'EMP4SDQ' in line:
                    mp4sdqenergy = float(line.split()[2])
                    mpenergies.append(mp4sdqenergy)
                # This is really MP4SD(T)Q.
                if 'EMP4 ' in line:
                    mp4sdtqenergy = float(line.split()[2])
                    mpenergies.append(mp4sdtqenergy)

                line = next(inputfile)

            mpenergies = [utils.convertor(mpe, 'hartree', 'eV')
                          for mpe in mpenergies]
            self.mpenergies.append(mpenergies)

        # Coupled cluster corrections.
        # Hopefully we only have to deal with ccman2 here.

        if 'CCD total energy' in line:
            if not hasattr(self, 'ccenergies'):
                self.ccenergies = []
            ccdenergy = float(line.split()[-1])
            ccdenergy = utils.convertor(ccdenergy, 'hartree', 'eV')
            self.ccenergies.append(ccdenergy)
        if 'CCSD total energy' in line:
            has_triples = False
            if not hasattr(self, 'ccenergies'):
                self.ccenergies = []
            ccsdenergy = float(line.split()[-1])
            # Make sure we aren't actually doing CCSD(T).
            line = next(inputfile)
            line = next(inputfile)
            if 'CCSD(T) total energy' in line:
                has_triples = True
                ccsdtenergy = float(line.split()[-1])
                ccsdtenergy = utils.convertor(ccsdtenergy, 'hartree', 'eV')
                self.ccenergies.append(ccsdtenergy)
            if not has_triples:
                ccsdenergy = utils.convertor(ccsdenergy, 'hartree', 'eV')
                self.ccenergies.append(ccsdenergy)

        # Electronic transitions. Works for both CIS and TDDFT.

        if 'Excitation Energies' in line:

            # Restricted:
            # ---------------------------------------------------
            #         TDDFT/TDA Excitation Energies
            # ---------------------------------------------------
            #
            # Excited state   1: excitation energy (eV) =    3.6052
            #    Total energy for state   1:   -382.167872200685
            #    Multiplicity: Triplet
            #    Trans. Mom.:  0.0000 X   0.0000 Y   0.0000 Z
            #    Strength   :  0.0000
            #    D( 33) --> V(  3) amplitude =  0.2618
            #    D( 34) --> V(  2) amplitude =  0.2125
            #    D( 35) --> V(  1) amplitude =  0.9266
            #
            # Unrestricted:
            # Excited state   2: excitation energy (eV) =    2.3156
            #    Total energy for state   2:   -381.980177630969
            #    <S**2>     :  0.7674
            #    Trans. Mom.: -2.7680 X  -0.1089 Y   0.0000 Z
            #    Strength   :  0.4353
            #    S(  1) --> V(  1) amplitude = -0.3105 alpha
            #    D( 34) --> S(  1) amplitude =  0.9322 beta

            self.skip_lines(inputfile, ['dashes', 'blank'])
            line = next(inputfile)

            etenergies = []
            etsyms = []
            etoscs = []
            etsecs = []
            spinmap = {'alpha': 0, 'beta': 1}

            while list(set(line.strip())) != ['-']:

                # Take the total energy for the state and subtract from the
                # ground state energy, rather than just the EE;
                # this will be more accurate.
                if 'Total energy for state' in line:
                    energy = utils.convertor(float(line.split()[-1]), 'hartree', 'cm-1')
                    etenergy = energy - utils.convertor(self.scfenergies[-1], 'eV', 'cm-1')
                    etenergies.append(etenergy)
                # if 'excitation energy' in line:
                #     etenergy = utils.convertor(float(line.split()[-1]), 'eV', 'cm-1')
                #     etenergies.append(etenergy)
                if 'Multiplicity' in line:
                    etsym = line.split()[1]
                    etsyms.append(etsym)
                if 'Strength' in line:
                    strength = float(line.split()[-1])
                    etoscs.append(strength)
                # This is the list of transitions.
                if 'amplitude' in line:
                    sec = []
                    while line.strip() != '':
                        if self.unrestricted:
                            spin = spinmap[line[42:47].strip()]
                        else:
                            spin = 0
                        startidx = int(line[6:9]) - 1
                        start = (startidx, spin)
                        # Q-Chem starts reindexing virtual orbitals at 1.
                        endidx = int(line[17:20]) - 1 + self.nalpha
                        end = (endidx, spin)
                        contrib = float(line[34:41].strip())
                        sec.append([start, end, contrib])
                        line = next(inputfile)
                    etsecs.append(sec)

                line = next(inputfile)

            if not hasattr(self, 'etenergies'):
                self.etenergies = numpy.array(etenergies)
            if not hasattr(self, 'etsyms'):
                self.etsyms = etsyms
            if not hasattr(self, 'etosecs'):
                self.etoscs = numpy.array(etoscs)
            if not hasattr(self, 'etsecs') and len(etsecs) > 0:
                self.etsecs = etsecs

        # Molecular orbital energies and symmetries.

        if 'Orbital Energies (a.u.) and Symmetries' in line:

            #  --------------------------------------------------------------
            #              Orbital Energies (a.u.) and Symmetries
            #  --------------------------------------------------------------
            #
            #  Alpha MOs, Restricted
            #  -- Occupied --
            # -10.018 -10.018 -10.008 -10.008 -10.007 -10.007 -10.006 -10.005
            #   1 Bu    1 Ag    2 Bu    2 Ag    3 Bu    3 Ag    4 Bu    4 Ag
            #  -9.992  -9.992  -0.818  -0.755  -0.721  -0.704  -0.670  -0.585
            #   5 Ag    5 Bu    6 Ag    6 Bu    7 Ag    7 Bu    8 Bu    8 Ag
            #  -0.561  -0.532  -0.512  -0.462  -0.439  -0.410  -0.400  -0.397
            #   9 Ag    9 Bu   10 Ag   11 Ag   10 Bu   11 Bu   12 Bu   12 Ag
            #  -0.376  -0.358  -0.349  -0.330  -0.305  -0.295  -0.281  -0.263
            #  13 Bu   14 Bu   13 Ag    1 Au   15 Bu   14 Ag   15 Ag    1 Bg
            #  -0.216  -0.198  -0.160
            #   2 Au    2 Bg    3 Bg
            #  -- Virtual --
            #   0.050   0.091   0.116   0.181   0.280   0.319   0.330   0.365
            #   3 Au    4 Au    4 Bg    5 Au    5 Bg   16 Ag   16 Bu   17 Bu
            #   0.370   0.413   0.416   0.422   0.446   0.469   0.496   0.539
            #  17 Ag   18 Bu   18 Ag   19 Bu   19 Ag   20 Bu   20 Ag   21 Ag
            #   0.571   0.587   0.610   0.627   0.646   0.693   0.743   0.806
            #  21 Bu   22 Ag   22 Bu   23 Bu   23 Ag   24 Ag   24 Bu   25 Ag
            #   0.816
            #  25 Bu
            #
            #  Beta MOs, Restricted
            #  -- Occupied --
            # -10.018 -10.018 -10.008 -10.008 -10.007 -10.007 -10.006 -10.005
            #   1 Bu    1 Ag    2 Bu    2 Ag    3 Bu    3 Ag    4 Bu    4 Ag
            #  -9.992  -9.992  -0.818  -0.755  -0.721  -0.704  -0.670  -0.585
            #   5 Ag    5 Bu    6 Ag    6 Bu    7 Ag    7 Bu    8 Bu    8 Ag
            #  -0.561  -0.532  -0.512  -0.462  -0.439  -0.410  -0.400  -0.397
            #   9 Ag    9 Bu   10 Ag   11 Ag   10 Bu   11 Bu   12 Bu   12 Ag
            #  -0.376  -0.358  -0.349  -0.330  -0.305  -0.295  -0.281  -0.263
            #  13 Bu   14 Bu   13 Ag    1 Au   15 Bu   14 Ag   15 Ag    1 Bg
            #  -0.216  -0.198  -0.160
            #   2 Au    2 Bg    3 Bg
            #  -- Virtual --
            #   0.050   0.091   0.116   0.181   0.280   0.319   0.330   0.365
            #   3 Au    4 Au    4 Bg    5 Au    5 Bg   16 Ag   16 Bu   17 Bu
            #   0.370   0.413   0.416   0.422   0.446   0.469   0.496   0.539
            #  17 Ag   18 Bu   18 Ag   19 Bu   19 Ag   20 Bu   20 Ag   21 Ag
            #   0.571   0.587   0.610   0.627   0.646   0.693   0.743   0.806
            #  21 Bu   22 Ag   22 Bu   23 Bu   23 Ag   24 Ag   24 Bu   25 Ag
            #   0.816
            #  25 Bu
            #  --------------------------------------------------------------

            self.skip_line(inputfile, 'dashes')
            line = next(inputfile)
            # Sometimes Q-Chem gets a little confused...
            while 'Warning : Irrep of orbital' in line:
                line = next(inputfile)
            line = next(inputfile)
            energies_alpha = []
            symbols_alpha = []
            if self.unrestricted:
                energies_beta = []
                symbols_beta = []
            line = next(inputfile)

            # The end of the block is either a blank line or only dashes.
            while not self.dashes_and_spaces.search(line):
                if 'Occupied' in line or 'Virtual' in line:
                    # A nice trick to find where the HOMO is.
                    if 'Virtual' in line:
                        self.homos = [len(energies_alpha)-1]
                    line = next(inputfile)
                # Parse the energies and symmetries in pairs of lines.
                # energies = [utils.convertor(energy, 'hartree', 'eV')
                #             for energy in map(float, line.split())]
                # This convoluted bit handles '*******' when present.
                energies = []
                energy_line = line.split()
                for e in energy_line:
                    try:
                        energy = utils.convertor(self.float(e), 'hartree', 'eV')
                    except ValueError:
                        energy = numpy.nan
                    energies.append(energy)
                energies_alpha.extend(energies)
                line = next(inputfile)
                symbols = line.split()[1::2]
                symbols_alpha.extend(symbols)
                line = next(inputfile)

            # Only look at the second block if doing an unrestricted calculation.
            # This might be a problem for ROHF/ROKS.
            if self.unrestricted:
                self.skip_line(inputfile, 'header')
                line = next(inputfile)
                while not self.dashes_and_spaces.search(line):
                    if 'Occupied' in line or 'Virtual' in line:
                        # This will definitely exist, thanks to the above block.
                        if 'Virtual' in line:
                            if len(self.homos) == 1:
                                self.homos.append(len(energies_beta)-1)
                        line = next(inputfile)
                    energies = []
                    energy_line = line.split()
                    for e in energy_line:
                        try:
                            energy = utils.convertor(self.float(e), 'hartree', 'eV')
                        except ValueError:
                            energy = numpy.nan
                        energies.append(energy)
                    energies_beta.extend(energies)
                    line = next(inputfile)
                    symbols = line.split()[1::2]
                    symbols_beta.extend(symbols)
                    line = next(inputfile)

            # For now, only keep the last set of MO energies, even though it is
            # printed at every step of geometry optimizations and fragment jobs.
            self.moenergies = [[]]
            self.mosyms = [[]]
            self.moenergies[0] = numpy.array(energies_alpha)
            self.mosyms[0] = symbols_alpha
            if self.unrestricted:
                self.moenergies.append([])
                self.mosyms.append([])
                self.moenergies[1] = numpy.array(energies_beta)
                self.mosyms[1] = symbols_beta
            self.set_attribute('nmo', len(self.moenergies[0]))

        # Molecular orbital energies, no symmetries.

        if line.strip() == 'Orbital Energies (a.u.)':

            # In the case of no orbital symmetries, the beta spin block is not
            # present for restricted calculations.

            #  --------------------------------------------------------------
            #                     Orbital Energies (a.u.)
            #  --------------------------------------------------------------
            #
            #  Alpha MOs
            #  -- Occupied --
            # ******* -38.595 -34.580 -34.579 -34.578 -19.372 -19.372 -19.364
            # -19.363 -19.362 -19.362  -4.738  -3.252  -3.250  -3.250  -1.379
            #  -1.371  -1.369  -1.365  -1.364  -1.362  -0.859  -0.855  -0.849
            #  -0.846  -0.840  -0.836  -0.810  -0.759  -0.732  -0.729  -0.704
            #  -0.701  -0.621  -0.610  -0.595  -0.587  -0.584  -0.578  -0.411
            #  -0.403  -0.355  -0.354  -0.352
            #  -- Virtual --
            #  -0.201  -0.117  -0.099  -0.086   0.020   0.031   0.055   0.067
            #   0.075   0.082   0.086   0.092   0.096   0.105   0.114   0.148
            #
            #  Beta MOs
            #  -- Occupied --
            # ******* -38.561 -34.550 -34.549 -34.549 -19.375 -19.375 -19.367
            # -19.367 -19.365 -19.365  -4.605  -3.105  -3.103  -3.102  -1.385
            #  -1.376  -1.376  -1.371  -1.370  -1.368  -0.863  -0.858  -0.853
            #  -0.849  -0.843  -0.839  -0.818  -0.765  -0.738  -0.737  -0.706
            #  -0.702  -0.624  -0.613  -0.600  -0.591  -0.588  -0.585  -0.291
            #  -0.291  -0.288  -0.275
            #  -- Virtual --
            #  -0.139  -0.122  -0.103   0.003   0.014   0.049   0.049   0.059
            #   0.061   0.070   0.076   0.081   0.086   0.090   0.098   0.106
            #   0.138
            #  --------------------------------------------------------------

            self.skip_lines(inputfile, ['dashes', 'blank'])
            line = next(inputfile)
            energies_alpha = []
            if self.unrestricted:
                energies_beta = []
            line = next(inputfile)

            # The end of the block is either a blank line or only dashes.
            while not self.dashes_and_spaces.search(line):
                if 'Occupied' in line or 'Virtual' in line:
                    # A nice trick to find where the HOMO is.
                    if 'Virtual' in line:
                        self.homos = [len(energies_alpha)-1]
                    line = next(inputfile)
                energies = []
                energy_line = line.split()
                for e in energy_line:
                    try:
                        energy = utils.convertor(self.float(e), 'hartree', 'eV')
                    except ValueError:
                        energy = numpy.nan
                    energies.append(energy)
                energies_alpha.extend(energies)
                line = next(inputfile)

            line = next(inputfile)
            # Only look at the second block if doing an unrestricted calculation.
            # This might be a problem for ROHF/ROKS.
            if self.unrestricted:
                self.skip_lines(inputfile, ['blank'])
                line = next(inputfile)
                while not self.dashes_and_spaces.search(line):
                    if 'Occupied' in line or 'Virtual' in line:
                        # This will definitely exist, thanks to the above block.
                        if 'Virtual' in line:
                            if len(self.homos) == 1:
                                self.homos.append(len(energies_beta)-1)
                        line = next(inputfile)
                    energies = []
                    energy_line = line.split()
                    for e in energy_line:
                        try:
                            energy = utils.convertor(self.float(e), 'hartree', 'eV')
                        except ValueError:
                            energy = numpy.nan
                        energies.append(energy)
                    energies_beta.extend(energies)
                    line = next(inputfile)

            # For now, only keep the last set of MO energies, even though it is
            # printed at every step of geometry optimizations and fragment jobs.
            self.moenergies = [[]]
            self.moenergies[0] = numpy.array(energies_alpha)
            if self.unrestricted:
                self.moenergies.append([])
                self.moenergies[1] = numpy.array(energies_beta)
            self.set_attribute('nmo', len(self.moenergies[0]))

        # Population analysis.

        if 'Ground-State Mulliken Net Atomic Charges' in line:
            self.parse_charge_section(inputfile, 'mulliken')
        if 'Hirshfeld Atomic Charges' in line:
            self.parse_charge_section(inputfile, 'hirshfeld')
        if 'Ground-State ChElPG Net Atomic Charges' in line:
            self.parse_charge_section(inputfile, 'chelpg')

        # Multipole moments are not printed in lexicographical order,
        # so we need to parse and sort them. The units seem OK, but there
        # is some uncertainty about the reference point and whether it
        # can be changed.
        #
        # Notice how the letter/coordinate labels change to coordinate ranks
        # after hexadecapole moments, and need to be translated. Additionally,
        # after 9-th order moments the ranks are not necessarily single digits
        # and so there are spaces between them.
        #
        # -----------------------------------------------------------------
        #                    Cartesian Multipole Moments
        #                  LMN = < X^L Y^M Z^N >
        # -----------------------------------------------------------------
        #    Charge (ESU x 10^10)
        #                 0.0000
        #    Dipole Moment (Debye)
        #         X       0.0000      Y       0.0000      Z       0.0000
        #       Tot       0.0000
        #    Quadrupole Moments (Debye-Ang)
        #        XX     -50.9647     XY      -0.1100     YY     -50.1441
        #        XZ       0.0000     YZ       0.0000     ZZ     -58.5742
        # ...
        #    5th-Order Moments (Debye-Ang^4)
        #       500       0.0159    410      -0.0010    320       0.0005
        #       230       0.0000    140       0.0005    050       0.0012
        # ...
        # -----------------------------------------------------------------
        #
        if "Cartesian Multipole Moments" in line:

            # This line appears not by default, but only when
            # `multipole_order` > 4:
            line = inputfile.next()
            if 'LMN = < X^L Y^M Z^N >' in line:
                line = inputfile.next()

            # The reference point is always the origin, although normally the molecule
            # is moved so that the center of charge is at the origin.
            self.reference = [0.0, 0.0, 0.0]
            self.moments = [self.reference]

            # Watch out! This charge is in statcoulombs without the exponent!
            # We should expect very good agreement, however Q-Chem prints
            # the charge only with 5 digits, so expect 1e-4 accuracy.
            charge_header = inputfile.next()
            assert charge_header.split()[0] == "Charge"
            charge = float(inputfile.next().strip())
            charge = utils.convertor(charge, 'statcoulomb', 'e') * 1e-10
            # Allow this to change until fragment jobs are properly implemented.
            # assert abs(charge - self.charge) < 1e-4

            # This will make sure Debyes are used (not sure if it can be changed).
            line = inputfile.next()
            assert line.strip() == "Dipole Moment (Debye)"

            while "-----" not in line:

                # The current multipole element will be gathered here.
                multipole = []

                line = inputfile.next()
                while ("-----" not in line) and ("Moment" not in line):

                    cols = line.split()

                    # The total (norm) is printed for dipole but not other multipoles.
                    if cols[0] == 'Tot':
                        line = inputfile.next()
                        continue

                    # Find and replace any 'stars' with NaN before moving on.
                    for i in range(len(cols)):
                        if '***' in cols[i]:
                            cols[i] = numpy.nan

                    # The moments come in pairs (label followed by value) up to the 9-th order,
                    # although above hexadecapoles the labels are digits representing the rank
                    # in each coordinate. Above the 9-th order, ranks are not always single digits,
                    # so there are spaces between them, which means moments come in quartets.
                    if len(self.moments) < 5:
                        for i in range(len(cols)//2):
                            lbl = cols[2*i]
                            m = cols[2*i + 1]
                            multipole.append([lbl, m])
                    elif len(self.moments) < 10:
                        for i in range(len(cols)//2):
                            lbl = cols[2*i]
                            lbl = 'X'*int(lbl[0]) + 'Y'*int(lbl[1]) + 'Z'*int(lbl[2])
                            m = cols[2*i + 1]
                            multipole.append([lbl, m])
                    else:
                        for i in range(len(cols)//4):
                            lbl = 'X'*int(cols[4*i]) + 'Y'*int(cols[4*i + 1]) + 'Z'*int(cols[4*i + 2])
                            m = cols[4*i + 3]
                            multipole.append([lbl, m])

                    line = inputfile.next()

                # Sort should use the first element when sorting lists,
                # so this should simply work, and afterwards we just need
                # to extract the second element in each list (the actual moment).
                multipole.sort()
                multipole = [m[1] for m in multipole]
                self.moments.append(multipole)

        # For `method = force` or geometry optimizations,
        # the gradient is printed.
        if 'Gradient of SCF Energy' in line:
            if not hasattr(self, 'grads'):
                self.grads = []
            grad = numpy.empty(shape=(3, self.natom))
            # A maximum of 6 columns/block.
            ncols = 6
            line = next(inputfile)
            colcounter = 0
            while colcounter < self.natom:
                if line[:5].strip() == '':
                    line = next(inputfile)
                rowcounter = 0
                while rowcounter < 3:
                    row = list(map(float, line.split()[1:]))
                    grad[rowcounter][colcounter:colcounter+ncols] = row
                    line = next(inputfile)
                    rowcounter += 1
                colcounter += ncols
            self.grads.append(grad.T)

        # For IR-related jobs, the Hessian is printed (dim: 3*natom, 3*natom).
        # Note that this is *not* the mass-weighted Hessian.
        if 'Hessian of the SCF Energy' in line:
            if not hasattr(self, 'hessian'):
                # A maximum of 6 columns/block.
                ncols = 6
                dim = 3*self.natom
                self.hessian = numpy.empty(shape=(dim, dim))
                line = next(inputfile)
                colcounter = 0
                while colcounter < dim:
                    if line[:5].strip() == '':
                        line = next(inputfile)
                    rowcounter = 0
                    while rowcounter < dim:
                        row = list(map(float, line.split()[1:]))
                        self.hessian[rowcounter][colcounter:colcounter+ncols] = row
                        line = next(inputfile)
                        rowcounter += 1
                    colcounter += ncols

        # Start of the IR/Raman frequency section.
        if 'VIBRATIONAL ANALYSIS' in line:

            while 'STANDARD THERMODYNAMIC QUANTITIES' not in line:

                ## IR, optional Raman:

                # **********************************************************************
                # **                                                                  **
                # **                       VIBRATIONAL ANALYSIS                       **
                # **                       --------------------                       **
                # **                                                                  **
                # **        VIBRATIONAL FREQUENCIES (CM**-1) AND NORMAL MODES         **
                # **     FORCE CONSTANTS (mDYN/ANGSTROM) AND REDUCED MASSES (AMU)     **
                # **                  INFRARED INTENSITIES (KM/MOL)                   **
                ##** RAMAN SCATTERING ACTIVITIES (A**4/AMU) AND DEPOLARIZATION RATIOS **
                # **                                                                  **
                # **********************************************************************


                # Mode:                 1                      2                      3
                # Frequency:      -106.88                -102.91                 161.77
                # Force Cnst:      0.0185                 0.0178                 0.0380
                # Red. Mass:       2.7502                 2.8542                 2.4660
                # IR Active:          NO                     YES                    YES
                # IR Intens:        0.000                  0.000                  0.419
                # Raman Active:       YES                    NO                     NO
                ##Raman Intens:     2.048                  0.000                  0.000
                ##Depolar:          0.750                  0.000                  0.000
                #               X      Y      Z        X      Y      Z        X      Y      Z
                # C          0.000  0.000 -0.100   -0.000  0.000 -0.070   -0.000 -0.000 -0.027
                # C          0.000  0.000  0.045   -0.000  0.000 -0.074    0.000 -0.000 -0.109
                # C          0.000  0.000  0.148   -0.000 -0.000 -0.074    0.000  0.000 -0.121
                # C          0.000  0.000  0.100   -0.000 -0.000 -0.070    0.000  0.000 -0.027
                # C          0.000  0.000 -0.045    0.000 -0.000 -0.074   -0.000 -0.000 -0.109
                # C          0.000  0.000 -0.148    0.000  0.000 -0.074   -0.000 -0.000 -0.121
                # H         -0.000  0.000  0.086   -0.000  0.000 -0.082    0.000 -0.000 -0.102
                # H          0.000  0.000  0.269   -0.000 -0.000 -0.091    0.000  0.000 -0.118
                # H          0.000  0.000 -0.086    0.000 -0.000 -0.082   -0.000  0.000 -0.102
                # H         -0.000  0.000 -0.269    0.000  0.000 -0.091   -0.000 -0.000 -0.118
                # C          0.000 -0.000  0.141   -0.000 -0.000 -0.062   -0.000  0.000  0.193
                # C         -0.000 -0.000 -0.160    0.000  0.000  0.254   -0.000  0.000  0.043
                # H          0.000 -0.000  0.378   -0.000  0.000 -0.289    0.000  0.000  0.519
                # H         -0.000 -0.000 -0.140    0.000  0.000  0.261   -0.000 -0.000  0.241
                # H         -0.000 -0.000 -0.422    0.000  0.000  0.499   -0.000  0.000 -0.285
                # C          0.000 -0.000 -0.141    0.000  0.000 -0.062   -0.000 -0.000  0.193
                # C         -0.000 -0.000  0.160   -0.000 -0.000  0.254    0.000  0.000  0.043
                # H          0.000 -0.000 -0.378    0.000 -0.000 -0.289   -0.000  0.000  0.519
                # H         -0.000 -0.000  0.140   -0.000 -0.000  0.261    0.000  0.000  0.241
                # H         -0.000 -0.000  0.422   -0.000 -0.000  0.499    0.000  0.000 -0.285
                # TransDip   0.000 -0.000 -0.000    0.000 -0.000 -0.000   -0.000  0.000  0.021

                # Mode:                 4                      5                      6
                # ...

                # There isn't any symmetry information for normal modes present
                # in Q-Chem.
                # if not hasattr(self, 'vibsyms'):
                #     self.vibsyms = []

                if 'Frequency:' in line:
                    if not hasattr(self, 'vibfreqs'):
                        self.vibfreqs = []
                    vibfreqs = map(float, line.split()[1:])
                    self.vibfreqs.extend(vibfreqs)

                if 'IR Intens:' in line:
                    if not hasattr(self, 'vibirs'):
                        self.vibirs = []
                    vibirs = map(float, line.split()[2:])
                    self.vibirs.extend(vibirs)

                if 'Raman Intens:' in line:
                    if not hasattr(self, 'vibramans'):
                        self.vibramans = []
                    vibramans = map(float, line.split()[2:])
                    self.vibramans.extend(vibramans)

                # This is the start of the displacement block.
                if line.split()[0:3] == ['X', 'Y', 'Z']:
                    if not hasattr(self, 'vibdisps'):
                        self.vibdisps = []
                    disps = []
                    for k in range(self.natom):
                        line = next(inputfile)
                        numbers = list(map(float, line.split()[1:]))
                        N = len(numbers) // 3
                        if not disps:
                            for n in range(N):
                                disps.append([])
                        for n in range(N):
                            disps[n].append(numbers[3*n:(3*n)+3])
                    self.vibdisps.extend(disps)

                line = next(inputfile)

                # Anharmonic vibrational analysis.
                # Q-Chem includes 3 theories: VPT2, TOSH, and VCI.
                # For now, just take the VPT2 results.

                # if 'VIBRATIONAL ANHARMONIC ANALYSIS' in line:

                #     while list(set(line.strip())) != ['=']:
                #         if 'VPT2' in line:
                #             if not hasattr(self, 'vibanharms'):
                #                 self.vibanharms = []
                #             self.vibanharms.append(float(line.split()[-1]))
                #         line = next(inputfile)

        if 'STANDARD THERMODYNAMIC QUANTITIES AT' in line:

            if not hasattr(self, 'temperature'):
                self.temperature = float(line.split()[4])
            # Not supported yet.
            if not hasattr(self, 'pressure'):
                self.pressure = float(line.split()[7])
            self.skip_lines(inputfile, ['blank', 'Imaginary'])
            line = next(inputfile)
            # Not supported yet.
            if 'Zero point vibrational energy' in line:
                if not hasattr(self, 'zpe'):
                    # Convert from kcal/mol to Hartree/particle.
                    self.zpe = utils.convertor(float(line.split()[4]),
                                               'kcal', 'hartree')

            atommasses = []

            while 'Archival summary' not in line:

                if 'Has Mass' in line:
                    atommass = float(line.split()[6])
                    atommasses.append(atommass)

                if 'Total Enthalpy' in line:
                    if not hasattr(self, 'enthalpy'):
                        enthalpy = float(line.split()[2])
                        self.enthalpy = utils.convertor(enthalpy,
                                                        'kcal', 'hartree')
                if 'Total Entropy' in line:
                    if not hasattr(self, 'entropy'):
                        entropy = float(line.split()[2]) * self.temperature / 1000
                        # This is the *temperature dependent* entropy.
                        self.entropy = utils.convertor(entropy,
                                                       'kcal', 'hartree')
                    if not hasattr(self, 'freeenergy'):
                        self.freeenergy = self.enthalpy - self.entropy

                line = next(inputfile)

            if not hasattr(self, 'atommasses'):
                self.atommasses = numpy.array(atommasses)

        # TODO:
        # 'aonames'
        # 'atombasis'
        # 'freeenergy'
        # 'fonames'
        # 'fooverlaps'
        # 'fragnames'
        # 'frags'
        # 'gbasis'
        # 'mocoeffs'
        # 'nocoeffs'
        # 'scancoords'
        # 'scanenergies'
        # 'scannames'
        # 'scanparm'

    def parse_charge_section(self, inputfile, chargetype):
        """Parse the population analysis charge block."""
        self.skip_line(inputfile, 'blank')
        line = next(inputfile)
        has_spins = False
        if 'Spin' in line:
            if not hasattr(self, 'atomspins'):
                self.atomspins = dict()
            has_spins = True
            spins = []
        self.skip_line(inputfile, 'dashes')
        if not hasattr(self, 'atomcharges'):
            self.atomcharges = dict()
        charges = []
        line = next(inputfile)

        while list(set(line.strip())) != ['-']:
            elements = line.split()
            charge = self.float(elements[2])
            charges.append(charge)
            if has_spins:
                spin = self.float(elements[3])
                spins.append(spin)
            line = next(inputfile)

        self.atomcharges[chargetype] = numpy.array(charges)
        if has_spins:
            self.atomspins[chargetype] = numpy.array(spins)


if __name__ == '__main__':
    import sys
    import doctest, qchemparser

    if len(sys.argv) == 1:
        doctest.testmod(qchemparser, verbose=False)

    if len(sys.argv) == 2:
        parser = qchemparser.QChem(sys.argv[1])
        data = parser.parse()

    if len(sys.argv) > 2:
        for i in range(len(sys.argv[2:])):
            if hasattr(data, sys.argv[2 + i]):
                print(getattr(data, sys.argv[2 + i]))
