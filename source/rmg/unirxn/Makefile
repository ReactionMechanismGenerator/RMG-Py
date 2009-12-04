################################################################################
#
#	Makefile for compiling rmg.unirxn Fortran modules
#
################################################################################

all: states.so mastereqn.so msc.so rs.so

states.so: states.f90
	f2py --fcompiler=$(F2PY_FCOMPILER) -c states.f90 -m states -lblas -llapack

mastereqn.so: mastereqn.f90
	f2py --fcompiler=$(F2PY_FCOMPILER) -c mastereqn.f90 -m mastereqn -lblas -llapack

msc.so: msc.f90
	f2py --fcompiler=$(F2PY_FCOMPILER) -c msc.f90 -lblas -llapack -m msc

rs.so: rs.f90
	f2py --fcompiler=$(F2PY_FCOMPILER) -c rs.f90 -lblas -llapack -m rs

clean:
	rm -rf *.so