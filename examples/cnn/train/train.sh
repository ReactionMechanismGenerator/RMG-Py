#!/bin/bash

# Make sure your PYTHONPATH includes the necessary paths to run RMG-Py 
# Uncomment the following lines if your PYTHONPATH is not stored in 
# your .bashrc file already. Be sure to modify the paths to the locations
# of your code as necessary.

# RMGPy=$HOME/code/RMG-Py/
# export PYTHONPATH=$PYTHONPATH:$RMGPy/
INPUT='predictor_input.py'
DATA_FILE='datasets.txt'
TRAIN_MODE='full_train'
BATCH_SIZE=1
source activate rmg_env
export KERAS_BACKEND=theano
python $RMGPy/scripts/train_cnn.py -i $INPUT -d ${DATA_FILE} -t ${TRAIN_MODE} -bs ${BATCH_SIZE}
source deactivate