# -*- coding: utf-8 -*-
#
# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2007-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Classes and tools for storing and handling parsed data"""


import numpy


class ccData(object):
    """Stores data extracted by cclib parsers

    Description of cclib attributes:
        aonames -- atomic orbital names (list of strings)
        aooverlaps -- atomic orbital overlap matrix (array[2])
        atombasis -- indices of atomic orbitals on each atom (list of lists)
        atomcharges -- atomic partial charges (dict of arrays[1])
        atomcoords -- atom coordinates (array[3], angstroms)
        atommasses -- atom masses (array[1], daltons)
        atomnos -- atomic numbers (array[1])
        atomspins -- atomic spin densities (dict of arrays[1])
        charge -- net charge of the system (integer)
        ccenergies -- molecular energies with Coupled-Cluster corrections (array[2], eV)
        coreelectrons -- number of core electrons in atom pseudopotentials (array[1])
        enthalpy -- sum of electronic and thermal enthalpies (float, hartree/particle)
        entropy -- entropy (float, hartree/particle)
        etenergies -- energies of electronic transitions (array[1], 1/cm)
        etoscs -- oscillator strengths of electronic transitions (array[1])
        etrotats -- rotatory strengths of electronic transitions (array[1], ??)
        etsecs -- singly-excited configurations for electronic transitions (list of lists)
        etsyms -- symmetries of electronic transitions (list of string)
        freeenergy -- sum of electronic and thermal free energies (float, hartree/particle)
        fonames -- fragment orbital names (list of strings)
        fooverlaps -- fragment orbital overlap matrix (array[2])
        fragnames -- names of fragments (list of strings)
        frags -- indices of atoms in a fragment (list of lists)
        gbasis -- coefficients and exponents of Gaussian basis functions (PyQuante format)
        geotargets -- targets for convergence of geometry optimization (array[1])
        geovalues -- current values for convergence of geometry optmization (array[1])
        grads -- current values of forces (gradients) in geometry optimization (array[3])
        hessian -- elements of the force constant matrix (array[1])
        homos -- molecular orbital indices of HOMO(s) (array[1])
        mocoeffs -- molecular orbital coefficients (list of arrays[2])
        moenergies -- molecular orbital energies (list of arrays[1], eV)
        moments -- molecular multipole moments (list of arrays[], a.u.)
        mosyms -- orbital symmetries (list of lists)
        mpenergies -- molecular electronic energies with Möller-Plesset corrections (array[2], eV)
        mult -- multiplicity of the system (integer)
        natom -- number of atoms (integer)
        nbasis -- number of basis functions (integer)
        nmo -- number of molecular orbitals (integer)
        nocoeffs -- natural orbital coefficients (array[2])
        optdone -- flags whether an optimization has converged (Boolean)
        scancoords -- geometries of each scan step (array[3], angstroms)
        scanenergies -- energies of potential energy surface (list)
        scannames -- names of varaibles scanned (list of strings)
        scanparm -- values of parameters in potential energy surface (list of tuples)
        scfenergies -- molecular electronic energies after SCF (Hartree-Fock, DFT) (array[1], eV)
        scftargets -- targets for convergence of the SCF (array[2])
        scfvalues -- current values for convergence of the SCF (list of arrays[2])
        temperature -- tempature used for Thermochemistry (float, kelvin)
        vibanharms -- vibrational anharmonicity constants (array[2], 1/cm)
        vibdisps -- cartesian displacement vectors (array[3], delta angstrom)
        vibfreqs -- vibrational frequencies (array[1], 1/cm)
        vibirs -- IR intensities (array[1], km/mol)
        vibramans -- Raman intensities (array[1], A^4/Da)
        vibsyms -- symmetries of vibrations (list of strings)
    (1) The term 'array' refers to a numpy array
    (2) The number of dimensions of an array is given in square brackets
    (3) Python indexes arrays/lists starting at zero, so if homos==[10], then
            the 11th molecular orbital is the HOMO
    """

    # The expected types for all supported attributes.
    _attrtypes = {
        "aonames":        list,
        "aooverlaps":     numpy.ndarray,
        "atombasis":      list,
        "atomcharges":    dict,
        "atomcoords":     numpy.ndarray,
        "atommasses":     numpy.ndarray,
        "atomnos":        numpy.ndarray,
        "atomspins":      dict,
        "ccenergies":     numpy.ndarray,
        "charge":         int,
        "coreelectrons":  numpy.ndarray,
        "enthalpy":       float,
        "entropy":        float,
        "etenergies":     numpy.ndarray,
        "etoscs":         numpy.ndarray,
        "etrotats":       numpy.ndarray,
        "etsecs":         list,
        "etsyms":         list,
        "freeenergy":     float,
        "fonames":        list,
        "fooverlaps":     numpy.ndarray,
        "fragnames":      list,
        "frags":          list,
        'gbasis':         list,
        "geotargets":     numpy.ndarray,
        "geovalues":      numpy.ndarray,
        "grads":          numpy.ndarray,
        "hessian":        numpy.ndarray,
        "homos":          numpy.ndarray,
        "mocoeffs":       list,
        "moenergies":     list,
        "molmass":        float,
        "moments":        list,
        "mosyms":         list,
        "mpenergies":     numpy.ndarray,
        "mult":           int,
        "natom":          int,
        "nbasis":         int,
        "nmo":            int,
        "nocoeffs":       numpy.ndarray,
        "optdone":        bool,
        "rotcons":        list,
        "rotsymm":	      int,
        "scancoords":     numpy.ndarray,
        "scanenergies":   list,
        "scannames":      list,
        "scanparm":       list,
        "scfenergies":    numpy.ndarray,
        "scftargets":     numpy.ndarray,
        "scfvalues":      list,
        "temperature":    float,
        "vibanharms":     numpy.ndarray,
        "vibdisps":       numpy.ndarray,
        "vibfreqs":       numpy.ndarray,
        "vibirs":         numpy.ndarray,
        "vibramans":      numpy.ndarray,
        "vibsyms":        list,
    }

    # The name of all attributes can be generated from the dictionary above.
    _attrlist = sorted(_attrtypes.keys())

    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ['atomnos', 'coreelectrons', 'homos']

    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'moments', 'scfvalues']
    
    # Attributes that should be dictionaries of arrays (double precision).
    _dictsofarrays = ["atomcharges", "atomspins"]

    def __init__(self, attributes={}):
        """Initialize the cclibData object.
        
        Normally called in the parse() method of a Logfile subclass.
        
        Inputs:
            attributes - optional dictionary of attributes to load as data
        """

        if attributes:
            self.setattributes(attributes)
        
    def listify(self):
        """Converts all attributes that are arrays or lists/dicts of arrays to lists."""
        
        attrlist = [k for k in self._attrlist if hasattr(self, k)]
        for k in attrlist:
            v = self._attrtypes[k]
            if v == numpy.ndarray:
                setattr(self, k, getattr(self, k).tolist())
            elif v == list and k in self._listsofarrays:
                setattr(self, k, [x.tolist() for x in getattr(self, k)])
            elif v == dict and k in self._dictsofarrays:
                items = getattr(self, k).iteritems()
                pairs = [(key, val.tolist()) for key, val in items]
                setattr(self, k, dict(pairs))
    
    def arrayify(self):
        """Converts appropriate attributes to arrays or lists/dicts of arrays."""
        
        attrlist = [k for k in self._attrlist if hasattr(self, k)]
        for k in attrlist:
            v = self._attrtypes[k]
            precision = 'd'
            if k in self._intarrays:
                precision = 'i'
            if v == numpy.ndarray:
                a = getattr(self, k)
                setattr(self, k, numpy.array(getattr(self, k), precision))
            elif v == list and k in self._listsofarrays:
                setattr(self, k, [numpy.array(x, precision) for x in getattr(self, k)])
            elif v == dict and k in self._dictsofarrays:
                items = getattr(self, k).items()
                pairs = [(key, numpy.array(val, precision)) for key, val in items]
                setattr(self, k, dict(pairs))

    def getattributes(self, tolists=False):
        """Returns a dictionary of existing data attributes.
        
        Inputs:
            tolists - flag to convert attributes to lists where applicable
        """
    
        if tolists:
            self.listify()
        attributes = {}
        for attr in self._attrlist:
            if hasattr(self, attr):
                attributes[attr] = getattr(self, attr)
        if tolists:
            self.arrayify()
        return attributes

    def setattributes(self, attributes):
        """Sets data attributes given in a dictionary.
        
        Inputs:
            attributes - dictionary of attributes to set
        Outputs:
            invalid - list of attributes names that were not set, which
                      means they are not specified in self._attrlist
        """
    
        if type(attributes) is not dict:
            raise TypeError("attributes must be in a dictionary")
    
        valid = [a for a in attributes if a in self._attrlist]
        invalid = [a for a in attributes if a not in self._attrlist]
    
        for attr in valid:
            setattr(self, attr, attributes[attr])

        self.arrayify()
        self.typecheck()

        return invalid

    def typecheck(self):
        """Check the types of all attributes.

        If an attribute does not match the expected type, then attempt to
        convert; if that fails, only then raise a TypeError.
        """

        self.arrayify()
        for attr in [a for a in self._attrlist if hasattr(self, a)]:

            val = getattr(self, attr)
            if type(val) == self._attrtypes[attr]:
                continue

            try:
                val = self._attrtypes[attr](val)
            except ValueError:
                args = (attr, type(val), self._attrtypes[attr])
                raise TypeError("attribute %s is %s instead of %s and could not be converted" % args)


class ccData_optdone_bool(ccData):
    """This is the version of ccData where optdone is a Boolean."""

    def __init__(self, *args, **kwargs):

        super(ccData_optdone_bool, self).__init__(*args, **kwargs)

        self._attrtypes['optdone'] = bool

    def setattributes(self, *args, **kwargs):

        invalid = super(ccData_optdone_bool, self).setattributes(*args, **kwargs)

        # Reduce optdone to a Boolean, because it will be parsed as a list. If this list has any element,
        # it means that there was an optimized structure and optdone should be True.
        if hasattr(self, 'optdone'):
            self.optdone = len(self.optdone) > 0
