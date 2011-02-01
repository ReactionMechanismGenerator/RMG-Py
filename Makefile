################################################################################
#
#   Makefile for RMG Py
#
################################################################################

all: 
	$(MAKE) -C chempy
	$(MAKE) -C measure
	$(MAKE) -C statesfit
	$(MAKE) -C solver

clean:
	$(MAKE) -C chempy clean
	$(MAKE) -C measure clean
	$(MAKE) -C statesfit clean
	$(MAKE) -C solver clean

cleanall:
	$(MAKE) -C chempy cleanall
	$(MAKE) -C measure cleanall
	$(MAKE) -C statesfit cleanall
	$(MAKE) -C solver cleanall
