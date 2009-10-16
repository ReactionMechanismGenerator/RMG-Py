################################################################################
#
#	Makefile to build Cython modules into compiled libraries
#
################################################################################

all: rmg/graph.so rmg/chem.so
	python setup.py build_ext --inplace

rmg/graph.so: rmg/graph.c

rmg/graph.c: rmg/graph.py rmg/graph.pxd rmg/chem.pxd
	cython rmg/graph.py

rmg/chem.so: rmg/chem.c

rmg/chem.c: rmg/chem.py rmg/chem.pxd
	cython rmg/chem.py
	
clean:
	rm -rf rmg/*.so rmg/*.c build/