#!/bin/bash

# Make sure your PYTHONPATH includes the necessary paths to run RMG-Py 
# Uncomment the following lines if your PYTHONPATH is not stored in 
# your .bashrc file already. Be sure to modify the paths to the locations
# of your code as necessary.

# export PYTHONPATH=$PYTHONPATH:$HOME/RMG-Py/
# export PYTHONPATH=$PYTHONPATH:$HOME/RMG-database/
# export PYTHONPATH=$PYTHONPATH:$HOME/PyDAS/
# export PYTHONPATH=$PYTHONPATH:$HOME/PyDQED/

# In order for the solvers to work properly, be sure to compile the PyDAS and PyDQED codes
# prior to running RMG.  Also be sure to compile RMG.  Sometimes with large code changes,
# you must first 'make clean' to removed old compiled files and recompile with the 'make'
# command.

# Run RMG on the input.py file.
source activate rmg_env
export KERAS_BACKEND=theano
python ../../../rmg.py input.py
source deactivate
