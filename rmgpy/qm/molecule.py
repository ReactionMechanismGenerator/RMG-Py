#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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

import logging
import math
import os

import numpy as np
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    logging.debug("To use QM calculations you must correctly install rdkit.")


import rmgpy.molecule
import rmgpy.qm.qmdata as qmdata
import rmgpy.qm.symmetry as symmetry
import rmgpy.quantity
import rmgpy.statmech
import rmgpy.thermo
from rmgpy.qm.qmdata import parse_cclib_data
from rmgpy.thermo import ThermoData


class Geometry(object):
    """
    A geometry, used for quantum calculations.
    
    Created from a molecule. Geometry estimated by RDKit.
    
    The attributes are:
    
    =================== ======================= ====================================
    Attribute           Type                    Description
    =================== ======================= ====================================
    `settings`          :class:`QMSettings`     Settings for QM calculations
    `unique_id`          ``str``                 A short ID such as an augmented InChI Key
    `molecule`          :class:`Molecule`       RMG Molecule object
    `unique_id_long`      ``str``                 A long, truly unique ID such as an augmented InChI
    =================== ======================= ====================================
    
    """

    def __init__(self, settings, unique_id, molecule, unique_id_long=None):
        self.settings = settings
        #: A short unique ID such as an augmented InChI Key.
        self.unique_id = unique_id
        self.molecule = molecule
        if unique_id_long is None:
            self.unique_id_long = unique_id
        else:
            #: Long, truly unique, ID, such as the augmented InChI.
            self.unique_id_long = unique_id_long

        # ToDo: why do we copy self.settings.fileStore into self.fileStore ?
        # (and same for .scratchDirectory)
        if self.settings:
            self.fileStore = self.settings.fileStore
            self.scratchDirectory = self.settings.scratchDirectory
        else:
            self.fileStore = None
            self.scratchDirectory = None

        if self.fileStore and not os.path.exists(self.fileStore):
            logging.info("Creating permanent directory %s for qm files." % os.path.abspath(self.fileStore))
            os.makedirs(self.fileStore)

        if self.scratchDirectory and not os.path.exists(self.scratchDirectory):
            logging.info("Creating scratch directory %s for qm files." % os.path.abspath(self.scratchDirectory))
            os.makedirs(self.scratchDirectory)

    def get_file_path(self, extension, scratch=True):
        """
        Returns the path to the file with the given extension.
        
        The provided extension should include the leading dot.
        If called with `scratch=False` then it will be in the `fileStore` directory,
        else `scratch=True` is assumed and it will be in the `scratchDirectory` directory.
        """
        return os.path.join(
            self.settings.scratchDirectory if scratch else self.settings.fileStore,
            self.unique_id + extension
        )

    def get_crude_mol_file_path(self):
        """Returns the path of the crude mol file."""
        return self.get_file_path('.crude.mol')

    def get_refined_mol_file_path(self):
        """Returns the path the the refined mol file."""
        return self.get_file_path('.refined.mol')

    def generate_rdkit_geometries(self):
        """
        Use RDKit to guess geometry.

        Save mol files of both crude and refined.
        Saves coordinates on atoms.
        """
        rdmol, rdatom_idx = self.rd_build()

        atoms = len(self.molecule.atoms)
        dist_geom_attempts = 1
        if atoms > 3:  # this check prevents the number of attempts from being negative
            dist_geom_attempts = 5 * (atoms - 3)  # number of conformer attempts is just a linear scaling with molecule size, due to time considerations in practice, it is probably more like 3^(n-3) or something like that

        rdmol, min_e_id = self.rd_embed(rdmol, dist_geom_attempts)
        self.save_coordinates_from_rdmol(rdmol, min_e_id, rdatom_idx)

    def rd_build(self):
        """
        Import rmg molecule and create rdkit molecule with the same atom labeling.
        """
        return self.molecule.to_rdkit_mol(remove_h=False, return_mapping=True)

    def rd_embed(self, rdmol, num_conf_attempts):
        """
        Embed the RDKit molecule and create the crude molecule file.
        """
        AllChem.EmbedMultipleConfs(rdmol, num_conf_attempts, randomSeed=1)

        energy = 0.0
        min_e_id = 0
        min_e = 9.999999e99  # start with a very high number, which would never be reached

        crude = Chem.Mol(rdmol.ToBinary())

        for i in range(rdmol.GetNumConformers()):
            AllChem.UFFOptimizeMolecule(rdmol, confId=i)
            energy = AllChem.UFFGetMoleculeForceField(rdmol, confId=i).CalcEnergy()
            if energy < min_e:
                min_e_id = i
                min_e = energy

        with open(self.get_crude_mol_file_path(), 'w') as out_3d_crude:
            out_3d_crude.write(Chem.MolToMolBlock(crude, confId=min_e_id))

        with open(self.get_refined_mol_file_path(), 'w') as out_3d:
            out_3d.write(Chem.MolToMolBlock(rdmol, confId=min_e_id))

        return rdmol, min_e_id

    def save_coordinates_from_rdmol(self, rdmol, min_e_id, rdatom_idx):
        # Save xyz coordinates on each atom in molecule ****
        for atom in self.molecule.atoms:
            point = rdmol.GetConformer(min_e_id).GetAtomPosition(atom.sorting_label)
            atom.coords = np.array([point.x, point.y, point.z])

    def save_coordinates_from_qm_data(self, qmdata):
        """
        Save geometry info from QMData (eg CCLibData)
        """
        raise NotImplementedError


def load_thermo_data_file(file_path):
    """
    Load the specified thermo data file and return the dictionary of its contents.
    
    Returns `None` if the file is invalid or missing.
    
    Checks that the returned dictionary contains at least InChI, adjacencyList, thermoData.
    """
    if not os.path.exists(file_path):
        return None
    try:
        with open(file_path) as result_file:
            logging.info('Reading existing thermo file {0}'.format(file_path))
            global_context = {'__builtins__': None}
            local_context = {
                '__builtins__': None,
                'True': True,
                'False': False,
                'ThermoData': rmgpy.thermo.ThermoData,
                'PointGroup': symmetry.PointGroup,
                'QMData': qmdata.QMData,
                'array': np.array,
                'int32': np.int32,
            }
            exec(result_file.read(), global_context, local_context)
    except IOError:
        logging.info("Couldn't read thermo file {0}".format(file_path))
        return None
    except (NameError, TypeError, SyntaxError) as e:
        logging.error('The thermo file "{0}" was invalid:'.format(file_path))
        logging.exception(e)
        return None
    if 'InChI' not in local_context:
        logging.error('The thermo file "{0}" did not contain an InChI.'.format(file_path))
        return None
    if 'adjacencyList' not in local_context:
        logging.error('The thermo file "{0}" did not contain adjacencyList.'.format(file_path))
        return None
    if 'thermoData' not in local_context:
        logging.error('The thermo file "{0}" did not contain thermoData.'.format(file_path))
        return None
    return local_context


class QMMolecule(object):
    """ 
    A base class for QM Molecule calculations.
    
    Specific programs and methods should inherit from this and define some
    extra attributes and methods:
    
     * outputFileExtension
     * inputFileExtension
     * generate_qm_data() ...and whatever else is needed to make this method work.
     
    The attributes are:
    
    =================== ======================= ====================================
    Attribute           Type                    Description
    =================== ======================= ====================================
    `molecule`          :class:`Molecule`       RMG Molecule object
    `settings`          :class:`QMSettings`     Settings for QM calculations
    `unique_id`         ``str``                 A short ID such as an augmented InChI Key
    `unique_id_long`    ``str``                 A long, truly unique ID such as an augmented InChI
    =================== ======================= ====================================
    
    """

    def __init__(self, molecule, settings):
        self.molecule = molecule
        self.settings = settings

        self.unique_id = self.molecule.to_augmented_inchi_key()
        self.unique_id_long = self.molecule.to_augmented_inchi()

    def get_file_path(self, extension, scratch=True):
        """
        Returns the path to the file with the given extension.
        
        The provided extension should include the leading dot.
        If called with `scratch=False` then it will be in the `fileStore` directory,
        else `scratch=True` is assumed and it will be in the `scratchDirectory` directory.
        """
        # ToDo: this is duplicated in Geometry class. Should be refactored.
        return os.path.join(
            self.settings.scratchDirectory if scratch else self.settings.fileStore,
            self.unique_id + extension
        )

    @property
    def output_file_path(self):
        """Get the output file name."""
        return self.get_file_path(self.outputFileExtension)

    @property
    def input_file_path(self):
        """Get the input file name."""
        return self.get_file_path(self.inputFileExtension)

    def get_thermo_file_path(self):
        """Returns the path the thermo data file."""
        return self.get_file_path('.thermo', scratch=False)

    @property
    def script_attempts(self):
        """The number of attempts with different script keywords"""
        return len(self.keywords)

    @property
    def max_attempts(self):
        """The total number of attempts to try"""
        return 2 * len(self.keywords)

    def initialize(self):
        """
        Do any startup tasks.
        """
        self.check_ready()

    def check_ready(self):
        """
        Check that it's ready to run calculations.
        """
        self.settings.check_all_set()
        self.check_paths()

    def check_paths(self):
        """
        Check the paths in the settings are OK. Make folders as necessary.
        """
        self.settings.fileStore = os.path.expandvars(self.settings.fileStore)  # to allow things like $HOME or $RMGpy
        self.settings.scratchDirectory = os.path.expandvars(self.settings.scratchDirectory)
        for path in [self.settings.fileStore, self.settings.scratchDirectory]:
            if not os.path.exists(path):
                logging.info("Creating directory %s for QM files." % os.path.abspath(path))
                os.makedirs(path)

    def create_geometry(self):
        """
        Creates self.geometry with RDKit geometries
        """
        self.geometry = Geometry(self.settings, self.unique_id, self.molecule, unique_id_long=self.unique_id_long)
        self.geometry.generate_rdkit_geometries()
        return self.geometry

    def parse(self):
        """
        Parses the results of the Mopac calculation, and returns a QMData object.
        """
        parser = self.get_parser(self.output_file_path)
        parser.logger.setLevel(logging.ERROR)  # cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
        parser.rotcons = [] # give it an attribute and it won't delete it, leaving it on the parser object
        parser.molmass = None # give it an attribute and it won't delete it, leaving it on the parser object
        cclib_data = parser.parse()
        radical_number = self.molecule.get_radical_count()
        cclib_data.rotcons = parser.rotcons # this hack required because rotcons not part of a default cclib data object
        cclib_data.molmass = parser.molmass # this hack required because rotcons not part of a default cclib data object
        qm_data = parse_cclib_data(cclib_data, radical_number + 1)  # Should `radical_number+1` be `self.molecule.multiplicity` in the next line of code? It's the electronic ground state degeneracy.
        return qm_data

    def generate_qm_data(self):
        """
        Calculate the QM data somehow and return a CCLibData object, or None if it fails.
        """
        raise NotImplementedError("This should be defined in a subclass that inherits from QMMolecule")
        return qmdata.QMData() or None

    def generate_thermo_data(self):
        """
        Generate Thermo Data via a QM calc. 
        
        Returns None if it fails.
        """
        self.initialize()

        # First, see if we already have it.
        if self.load_thermo_data():
            return self.thermo

        # If not, generate the QM data
        self.qm_data = self.generate_qm_data()

        # If that fails, give up and return None.
        if self.qm_data is None:
            return None

        self.determine_point_group()

        # If that fails, give up and return None.
        if self.point_group is None:
            return None

        self.calculate_thermo_data()
        Cp0 = self.molecule.calculate_cp0()
        CpInf = self.molecule.calculate_cpinf()
        self.thermo.Cp0 = (Cp0, "J/(mol*K)")
        self.thermo.CpInf = (CpInf, "J/(mol*K)")
        self.save_thermo_data()
        return self.thermo

    def save_thermo_data(self):
        """
        Save the generated thermo data.
        """
        self.thermo.H298.units = 'kcal/mol'
        self.thermo.S298.units = 'cal/mol/K'
        self.thermo.Cpdata.units = 'cal/mol/K'
        with open(self.get_thermo_file_path(), 'w') as result_file:
            result_file.write('InChI = "{0!s}"\n'.format(self.unique_id_long))
            result_file.write("thermoData = {0!r}\n".format(self.thermo))
            result_file.write("pointGroup = {0!r}\n".format(self.point_group))
            result_file.write("qmData = {0!r}\n".format(self.qm_data))
            result_file.write('adjacencyList = """\n{0!s}"""\n'.format(self.molecule.to_adjacency_list(remove_h=False)))

    def load_thermo_data(self):
        """
        Try loading a thermo data from a previous run.
        """
        file_path = self.get_thermo_file_path()
        local_context = load_thermo_data_file(file_path)
        if local_context is None:
            # file does not exist or is invalid
            return None
        if local_context['InChI'] != self.unique_id_long:
            logging.error('The InChI in the thermo file {0} did not match the '
                          'current molecule {1}'.format(file_path, self.unique_id_long))
            return None
        loaded_molecule = rmgpy.molecule.Molecule().from_adjacency_list(local_context['adjacencyList'])
        if not loaded_molecule.is_isomorphic(self.molecule):
            logging.error('The adjacencyList in thermo file {0} did not match the '
                          'current molecule {1}'.format(file_path, self.unique_id_long))
            return None
        thermo = local_context['thermoData']
        assert isinstance(thermo, rmgpy.thermo.ThermoData)
        self.thermo = thermo

        self.point_group = symmetry.point_group_dictionary[local_context['pointGroup'].point_group]  # point to the one in the module level dictionary with the same name
        self.qm_data = local_context['qmData']
        return thermo

    def get_augmented_inchi_key(self):
        """
        Returns the augmented InChI from self.molecule 
        """
        return self.molecule.to_augmented_inchi_key()

    def get_mol_file_path_for_calculation(self, attempt):
        """
        Get the path to the MOL file of the geometry to use for calculation `attempt`.
        
        If attempt <= self.script_attempts then we use the refined coordinates,
        then we start to use the crude coordinates.
        """
        if attempt <= self.script_attempts:  # use UFF-refined coordinates
            return self.geometry.get_refined_mol_file_path()
        else:
            return self.geometry.get_crude_mol_file_path()

    def determine_point_group(self):
        """
        Determine point group using the SYMMETRY Program
        
        Stores the resulting :class:`PointGroup` in self.point_group
        """
        assert self.qm_data, "Need QM Data first in order to calculate point group."
        pgc = symmetry.PointGroupCalculator(self.settings, self.unique_id, self.qm_data)
        self.point_group = pgc.calculate()

    def calculate_chirality_correction(self):
        """
        Returns the chirality correction to entropy (R*ln(2) if chiral) in J/mol/K.
        """
        if self.point_group.chiral:
            return rmgpy.quantity.constants.R * math.log(2)
        else:
            return 0.

    def calculate_thermo_data(self):
        """
        Calculate the thermodynamic properties.
        
        Stores and returns a ThermoData object as self.thermo.
        self.qm_data and self.point_group need to be generated before this method is called.
        """
        assert self.qm_data, "Need QM Data first in order to calculate thermo."
        assert self.point_group, "Need Point Group first in order to calculate thermo."

        mass = getattr(self.qm_data, 'molecularMass', None)
        if mass is None:
            # If using a cclib that doesn't read molecular mass, for example
            mass = sum(rmgpy.molecule.element.get_element(int(a)).mass for a in self.qm_data.atomicNumbers)
            mass = rmgpy.quantity.Mass(mass, 'kg/mol')
        trans = rmgpy.statmech.IdealGasTranslation(mass=mass)
        if self.point_group.linear:
            # there should only be one rotational constant for a linear rotor
            rotational_constant = rmgpy.quantity.Frequency(max(self.qm_data.rotationalConstants.value),
                                                           self.qm_data.rotationalConstants.units)
            rot = rmgpy.statmech.LinearRotor(
                rotationalConstant=rotational_constant,
                symmetry=self.point_group.symmetry_number,
            )
        else:
            rot = rmgpy.statmech.NonlinearRotor(
                rotationalConstant=self.qm_data.rotationalConstants,
                symmetry=self.point_group.symmetry_number,
            )
        # @todo: should we worry about spherical top rotors?
        vib = rmgpy.statmech.HarmonicOscillator(frequencies=self.qm_data.frequencies)

        # @todo: We need to extract or calculate E0 somehow from the qmdata
        E0 = (0, "kJ/mol")
        self.statesmodel = rmgpy.statmech.Conformer(E0=E0,
                                                    modes=[trans, rot, vib],
                                                    spin_multiplicity=self.qm_data.groundStateDegeneracy)

        # we will use number of atoms from above (alternatively, we could use the chemGraph); this is needed to test whether the species is monoatomic
        # SI units are J/mol, but converted to kJ/mol for generating the thermo.
        Hf298 = self.qm_data.energy.value_si / 1000

        S298 = self.statesmodel.get_entropy(298.0)
        Tdata = [300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0]
        Cp = [self.statesmodel.get_heat_capacity(T) for T in Tdata]
        S298 = S298 + self.calculate_chirality_correction()
        comment = self.qm_data.source or "QM calculation of some sort."

        thermo = ThermoData(
            Tdata=(Tdata, "K"),
            Cpdata=(Cp, "J/(mol*K)"),
            H298=(Hf298, "kJ/mol"),
            S298=(S298, "J/(mol*K)"),
            Tmin=(300.0, "K"),
            Tmax=(2000.0, "K"),
            comment=comment
        )
        self.thermo = thermo
        return thermo


def parse_gaussian_molmass_and_rotcons(outfile):
    """
    Temporary function to parse molecular mass and rotational constants from Gaussian files until cclib is updated

    Args:
        outfile (str): Location of Gaussian file to parse

    Returns:
        float, list: Molecular mass, list of rotational constants
    """
    molmass = None
    rotcons = []

    with open(outfile, 'r') as f:
        lines = f.readlines()

    for line in lines:
        # Extract Rotational Constants
        # Example:
        # Rotational constants (GHZ):           3.13081     1.24272     0.88960
        # OR for linear molecules:
        # Rotational constants (GHZ): ************ 12.73690 12.73690
        # Note: rotational constant is converted to wavenumber units (1/cm) to standardize across parsers
        if 'Rotational constants (GHZ)' in line:
            splits = line.split()

            # Determine if the molecule is linear and only has two constants
            if '*' in line:  # linear molecule
                rotcons.append([0.0] + [float(splits[i]) / 29.9792458 for i in (-2, -1)])
            else:
                rotcons.append([float(splits[i]) / 29.9792458 for i in (-3, -2, -1)])

        # Extract Molecular Mass (in amu)
        # Example:
        # Molecular mass:   128.06260 amu.
        if 'Molecular mass:' in line:
            splits = line.split()
            molmass = float(splits[2])

    return molmass, rotcons
