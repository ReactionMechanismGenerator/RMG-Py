################################################################################
#
#   Makefile for RMG Py
#
################################################################################

all: 
	$(MAKE) -C chempy
	$(MAKE) -C measure
	$(MAKE) -C statesfit

clean:
	$(MAKE) -C chempy clean
	$(MAKE) -C measure clean
	$(MAKE) -C statesfit clean

cleanall:
	$(MAKE) -C chempy cleanall
	$(MAKE) -C measure cleanall
	$(MAKE) -C statesfit cleanall
