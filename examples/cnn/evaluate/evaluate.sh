#!/bin/bash

# Make sure your PYTHONPATH includes the necessary paths to run RMG-Py 
# Uncomment the following lines if your PYTHONPATH is not stored in 
# your .bashrc file already. Be sure to modify the paths to the locations
# of your code as necessary.

# RMGPy=$HOME/code/RMG-Py/
# export PYTHONPATH=$PYTHONPATH:$RMGPy/

CNN_MODEL='test_model'
DATA_FILE='test_datasets.txt'
source activate rmg_env
python $RMGPy/evaluate_cnn.py -d ${DATA_FILE} -m ${CNN_MODEL}
source deactivate