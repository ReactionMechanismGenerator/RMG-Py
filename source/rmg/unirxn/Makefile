################################################################################
#
#	Makefile for compiling rmg.unirxn Fortran modules
#
################################################################################

all: states.so mastereqn.so msc.so rs.so

states.so: states.f90
	f2py --fcompiler=$(F2PY_FCOMPILER) -c states.f90 -m states $(LDFLAGS)

mastereqn.so: mastereqn.f90
	f2py --fcompiler=$(F2PY_FCOMPILER) -c mastereqn.f90 -m mastereqn $(LDFLAGS)

msc.so: msc.f90
	f2py --fcompiler=$(F2PY_FCOMPILER) -c msc.f90 -m msc $(LDFLAGS)

rs.so: rs.f90
	f2py --fcompiler=$(F2PY_FCOMPILER) -c rs.f90 -m rs $(LDFLAGS)

clean:
	rm -rf *.so