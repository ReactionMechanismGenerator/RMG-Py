################################################################################
#
#	Makefile to build Cython and Fortran modules into compiled libraries
#
################################################################################

#  To run with  GNU Fotrtan 95 (gfortran) you need to fix numpy distutils as follows:
#    Find the file python2.5/site-packages/numpy/distutils/fcompiler/gnu.py
#    (which may be somewhere obscure like /Library/Frameworks/Python.framework/Versions/2.5/lib/python2.5/site-packages/numpy-1.3.0-py2.5-macosx-10.3-fat.egg/numpy/distutils/fcompiler/gnu.py )
#    replace line line 227:
#          'linker_so'    : ["<F90>", "-Wall"],
#    with this:
#          'linker_so'    : ["<F90>", "-Wall -dynamiclib -L/opt/local/lib -lpython2.5"],
#F90=gfortran
#F2PY_FCOMPILER=gnu95
#F90FLAGS=-fPIC -r8 -fbounds-check -W -Wall -Werror # -ftrace=full 

# Richard currently can't get this to work, so just leaving Josh's g95 version hopefully unbroken:

# Alternatively, use G95:

MAKE=make

F90=g95
F2PY_FCOMPILER=g95
F90FLAGS=-fPIC -r8 -ftrace=full -fbounds-check -Wall

################################################################################

all: cython spectral unirxn
	python setup.py build_ext --inplace

cython: rmg/graph.py rmg/chem.py rmg/thermo.py 

spectral: rmg/spectral/fit.so

unirxn: rmg/unirxn/states.so rmg/unirxn/mastereqn.so rmg/unirxn/msc.so rmg/unirxn/rs.so

rmg/unirxn/%.so: rmg/unirxn/%.f90
	$(MAKE) -C rmg/unirxn F90=$(F90) F2PY_FCOMPILER=$(F2PY_FCOMPILER) F90FLAGS="$(F90FLAGS)"

rmg/spectral/%.so: rmg/spectral/%.f90
	$(MAKE) -C rmg/spectral F90=$(F90) F2PY_FCOMPILER=$(F2PY_FCOMPILER) F90FLAGS="$(F90FLAGS)" LDFLAGS=$(LDFLAGS)

clean:
	$(MAKE) -C rmg/spectral clean
	$(MAKE) -C rmg/unirxn clean
	rm -rf rmg/*.so rmg/*.c build/