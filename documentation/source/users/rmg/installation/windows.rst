.. _windows:

********************
Windows Installation
********************

.. NOTE::
	RMG-Py is currently only compatible with 32-bit systems

.. _path:
	
Path Adjustments
================

Before beginning the installations listed below, go ahead and adjust the PATH environment variable now. This can be found with the Run command: ::

	systempropertiesadvanced
	
Click on "Environment Variables," then edit the System Variable "Path" by appending the following to the variable value: ::

	C:\Program Files\Git\bin;C:\MinGW\bin;C:\Python27;C:\Python27\Scripts
	
.. _git:

Git
===

* Create a `Github account <https://github.com/signup/free>`_, then download and run the `Git Installer <http://git-scm.com/download/win>`_. Be sure to install to the default "C:\\Program Files\\Git" directory.

* Provide easy access to the Git Bash command shell that will be used throughout this installation by creating a shortcut with the following target: ::

	C:\Windows\System32\cmd.exe /c "sh --login -i"

if this target does not work try using this one
	C:\Windows\System32\cmd.exe /c  ""C:\Program Files\Git\bin\sh.exe" --login -i"

* Generate and view your SSH Key by typing this into the Git Bash command line: ::

	cd ~
	ssh-keygen		# press enter to save to the default directory
				# create a password if desired
	cat .ssh/id_rsa.pub

  Right-click the title of your bash window and select Edit > Mark. Highlight the entire block of text and press enter to copy the selection to your clipboard. Add this SSH Key to your Github account `here <https://github.com/settings/ssh>`_.

.. _mingw:

MinGW
=====

Download and run the `MinGW installer <http://hivelocity.dl.sourceforge.net/project/mingw/Installer/mingw-get-inst/mingw-get-inst-20120426/mingw-get-inst-20120426.exe>`_. When prompted, select the "C compiler" and "Fortran compiler" options. Be sure to install to the "C:\\MinGW" directory; installing to the "Program Files" directory may lead to build issues later on.

.. _python:

Python
======

* Download and run the `Python 2.7 installer <http://www.python.org/ftp/python/2.7.3/python-2.7.3.msi>`_. Be sure to install to the default "C:\\Python27" directory.

* Patch distutils by opening "C:\\Python27\\Lib\\disutils\\cygwinccompiler.py" in a text editor and commenting out lines 132-137 and 322-328 with a # symbol at the beginning of each line: ::

	[132:137]
	
	#        self.set_executables(compiler='gcc -mcygwin -O -Wall',
	#                             compiler_so='gcc -mcygwin -mdll -O -Wall',
	#                             compiler_cxx='g++ -mcygwin -O -Wall',
	#                             linker_exe='gcc -mcygwin',
	#                             linker_so=('%s -mcygwin %s' %
	#                                        (self.linker_dll, shared_option)))
	
	[322:328]
	
	#        self.set_executables(compiler='gcc -mno-cygwin -O -Wall',
	#                             compiler_so='gcc -mno-cygwin -mdll -O -Wall',
	#                             compiler_cxx='g++ -mno-cygwin -O -Wall',
	#                             linker_exe='gcc -mno-cygwin',
	#                             linker_so='%s -mno-cygwin %s %s'
	#                                        % (self.linker_dll, shared_option,
	#                                           entry_point))

* Point distutils to MinGW: ::

	echo -e "[build]\ncompiler=mingw32" > /c/Python27/Lib/distutils/distutils.cfg

* Add Python bindings to MinGW's library: ::

	mingw-get install pexports
	pexports $SYSTEMROOT/system32/python27.dll > python27.def
	dlltool -D python27.dll -d python27.def -l libpython27.a
	mv libpython27.a /c/MinGW/lib/libpython27.a
	rm python27.def

  If the ``pexports`` step doesn't work then you can download the :file:`python27.def` file from the link on the `Cython wiki <http://wiki.cython.org/InstallingOnWindows>`_ and continue from the ``dlltool`` step.

	
.. _prepackageddependencies:

Prepackaged Dependencies
========================

Download and run the installers listed below. These builds have been verified as compatible, and several of these binaries appear to be the only dependable way to build those dependencies.

* `NumPy 1.6.2 <http://softlayer.dl.sourceforge.net/project/numpy/NumPy/1.6.2/numpy-1.6.2-win32-superpack-python2.7.exe>`_
* `SciPy 0.10.1 <http://softlayer.dl.sourceforge.net/project/scipy/scipy/0.10.1/scipy-0.10.1-win32-superpack-python2.7.exe>`_
* `matplotlib  1.1.0 <http://softlayer.dl.sourceforge.net/project/matplotlib/matplotlib/matplotlib-1.1.0/matplotlib-1.1.0.win32-py2.7.exe>`_
* `guppy 0.1.10 <http://www.sistemasagiles.com.ar/soft/guppy-0.1.10.win32-py2.7.exe>`_
	If you have Norton Antivirus on your computer it may try to remove this after you install it
* `OpenBabel 2.3.1 <http://voxel.dl.sourceforge.net/project/openbabel/openbabel/2.3.1/OpenBabel2.3.1_Windows_Installer.exe>`_

  * The OpenBabel installer includes some libraries (.dll files) that you also need for other purposes, so copy them out of the OpenBabel program directory and into your system directory so they are generally accessible: ::
	
		cd /c/PROGRA~1/OpenBabel-2.3.1
		cp libcairo-2.dll libpng14-14.dll zlib1.dll $SYSTEMROOT/System32

* `openbabel-python 1.7 <http://softlayer.dl.sourceforge.net/project/openbabel/openbabel-python/1.7/openbabel-python-1.7.py27.exe>`_
* `py2cairo 1.10.0 <http://wxpython.org/cairo/py2cairo-1.10.0.win32-py2.7.exe>`_
* `Graphviz 2.28.0 <http://www.graphviz.org/pub/graphviz/stable/windows/graphviz-2.28.0.msi>`_

.. _remainingdependencies:


RDKit
================

Project home on GitHub: https://github.com/rdkit/rdkit

Installation instructions: http://code.google.com/p/rdkit/wiki/GettingStarted
Build it with InChI support.
Required environment variables:

* RDBASE pointing to the root of the distribution
* PYTHONPATH: should include $RDBASE
* PATH: should include $RDBASE/lib

Remaining Dependencies
======================

Install the remaining six python dependencies using 'pip': ::

	curl https://raw.github.com/pypa/pip/master/contrib/get-pip.py | python
	easy_install pip
	pip install nose quantities sphinx pydot psutil xlwt cython

.. _rmgsources:

RMG
===

* Download all RMG source packages: ::

	cd /c
	git clone git@github.com:jwallen/PyDAS.git
	git clone git@github.com:jwallen/PyDQED.git
	git clone git@github.com:GreenGroup/RMG-database.git
	git clone git@github.com:GreenGroup/RMG-Py.git

* Build PyDAS by running the provided "make.bat" file, then install it: ::

	cd /c/PyDAS
	python setup.py install

* Build and install PyDQED: ::

	cd /c/PyDQED
	mingw32-make
	python setup.py install

* Build and install RMG-Py: ::

	cd /c/RMG-Py
	mingw32-make

* Run an example: ::

	cd /c/RMG-Py
	python rmg.py examples/rmg/minimal/input.py

  Verify your installation by opening the resulting `output.html <file:///C:/RMG-Py/examples/rmg/minimal/output.html>`_ file.
