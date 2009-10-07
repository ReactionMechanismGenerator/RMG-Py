# Makefile to run cython on the necessary files

# If graph.so and chem.so up to date, nothing to do
all: rmg/graph.so rmg/chem.so
	python setup.py build_ext --inplace

# remember graph.so depens on chem.pxd
rmg/graph.so: rmg/graph.pyx rmg/graph.pxd rmg/chem.pxd
	rm -f rmg/rmg.c rmg/rmg.so # delete the intermediate files so cython rebuilds

rmg/chem.so: rmg/chem.pyx rmg/chem.pxd
	rm -f rmg/chem.c rmg/chem.so # delete the intermediate files so cython rebuilds
	
cgit:
	git checkout -- rmg/graph.c rmg/chem.c # kill local changes to c files
	