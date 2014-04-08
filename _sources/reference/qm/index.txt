*************************************************
QMTP (:mod:`rmgpy.qm`)
*************************************************

.. module:: rmgpy.qm

The :mod:`rmgpy.qm` subpackage contains classes and functions for working
with molecular geometries, and interfacing with quantum chemistry software.



Main
======

.. currentmodule:: rmgpy.qm.main

======================= ========================================================
Class                   Description
======================= ========================================================
:class:`QMSettings`     A class to store settings related to quantum mechanics calculations
:class:`QMCalculator`   An object to store settings and previous calculations
======================= ========================================================



Molecule
==========

.. currentmodule:: rmgpy.qm.molecule

======================= ========================================================
Class                   Description
======================= ========================================================
:class:`Geometry`       A geometry, used for quantum calculations
:class:`QMMolecule`     A base class for QM Molecule calculations
======================= ========================================================



QM Data
=========

.. currentmodule:: rmgpy.qm.qmdata

======================= ========================================================
Class/Function          Description
======================= ========================================================
:class:`QMData`         General class for data extracted from a QM calculation
:class:`CCLibData`      QM Data extracted from a cclib data object
======================= ========================================================


QM Verifier
=============

.. currentmodule:: rmgpy.qm.qmverifier

======================= ========================================================
Class/Function          Description
======================= ========================================================
:class:`QMVerifier`     Verifies whether a QM job was succesfully completed
======================= ========================================================


Symmetry
==========

.. currentmodule:: rmgpy.qm.symmetry

==============================	========================================================
Class/Function          			Description
==============================	========================================================
:class:`PointGroup`         		A symmetry Point Group
:class:`PointGroupCalculator`    Wrapper type to determine molecular symmetry point groups based on 3D coordinates
:class:`SymmetryJob`				Determine the point group using the SYMMETRY program 
==============================	========================================================


Gaussian
===========

.. currentmodule:: rmgpy.qm.gaussian

======================= ========================================================
Class/Function          Description
======================= ========================================================
:class:`Gaussian`       A base class for all QM calculations that use Gaussian
:class:`GaussianMol`    A base Class for calculations of molecules using Gaussian. 
:class:`GaussianMolPM3` A base Class for calculations of molecules using Gaussian at PM3. 
:class:`GaussianMolPM6` A base Class for calculations of molecules using Gaussian at PM6. 
======================= ========================================================


Mopac
=======

.. currentmodule:: rmgpy.qm.mopac

======================= ========================================================
Class/Function          Description
======================= ========================================================
:class:`Mopac`          A base class for all QM calculations that use Mopac
:class:`MopacMol`       A base Class for calculations of molecules using Mopac. 
:class:`MopacMolPM3`    A base Class for calculations of molecules using Mopac at PM3. 
:class:`MopacMolPM6`    A base Class for calculations of molecules using Mopac at PM6. 
:class:`MopacMolPM7`    A base Class for calculations of molecules using Mopac at PM7. 
======================= ========================================================



.. toctree::
	:hidden:

	main
	molecule
	qmdata
	qmverifier
	symmetry
	gaussian
	mopac
	