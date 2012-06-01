.. _linux:

******************
Linux Installation
******************

RMG-Py and all of its dependencies may be easily installed through a short series of Terminal commands. The instructions listed below have been confirmed on a fresh Ubuntu 12.04 installation and should generally apply to other distributions.

* Install compilers and libraries: ::

	sudo apt-get install git g++ gfortran python-dev liblapack-dev python-openbabel python-setuptools python-pip

* After creating a `Github account <https://github.com/signup/free>`_, generate your public key: ::

	cd ~; ssh-keygen		# press enter to save to the default directory
					# create a password if desired
	cat .ssh/id_rsa.pub

  Copy this public key to your `Github profile <https://github.com/settings/ssh>`_.

* Install dependencies: ::

	sudo pip install numpy		# install NumPy before other packages
	
	sudo pip install scipy cython nose matplotlib quantities guppy sphinx psutil xlwt
	
	cd ~
	git clone git@github.com:jwallen/PyDAS.git
	git clone git@github.com:jwallen/PyDQED.git
	cd PyDAS; make F77=gfortran; sudo make install
	cd ../PyDQED; make F77=gfortran; sudo make install

* Install RMG-Py: ::

	cd ~
	git clone git@github.com:GreenGroup/RMG-database.git
	git clone git@github.com:GreenGroup/RMG-Py.git
	cd RMG-Py; make

* Run an example: ::

	python rmg.py examples/rmg/minimal/input.py

  Verify your installation by opening the resulting output.html file under the "examples/rmg/minimal" directory.
