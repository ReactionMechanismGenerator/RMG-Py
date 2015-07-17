.. _QMthermoInstall:

*******************************
Quantum Mechanical Calculations
*******************************

To run the :ref:`QMTP interface <qm>`, or quantum mechanical calculations for improved thermochemistry estimates of cyclic species,
you will need some additional software.


RDKit
=====

Project home on GitHub: https://github.com/rdkit/rdkit

Installation instructions: http://code.google.com/p/rdkit/wiki/GettingStarted
Build it with InChI support.

Mac users with `homebrew <http://brew.sh/>`_ can install it most easily with::

	brew tap rdkit/rdkit
	brew install rdkit --with-inchi

You'll need various environment variables set, eg.::

	export RDBASE=$HOME/rdkit # CHECK THIS (maybe you put RDKit somewhere else)
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$RDBASE/lib
	export PYTHONPATH=$PYTHONPATH:$RDBASE


SYMMETRY
========

See http://www.cobalt.chem.ucalgary.ca/ps/symmetry.
The source code for this is included in the `external` folder, but it needs compiling.::

	$ make -C $RMGpy/external/symmetry

This is now included in the master Makefile, so you can simply type::

	$ make QM


MOPAC
=====

This is a semi-empirical software package you can use for QM calculations.
See http://openmopac.net/

Ensure your environment contains the variables `MOPAC_LICENSE` and `MOPAC_DIR`. eg.::

	export MOPAC_LICENSE=$HOME/mopac/
	export MOPAC_DIR=$HOME/mopac

Gaussian
========

This can be used instead of MOPAC for QM calculations.

See http://www.gaussian.com

