
from rmgpy.cnn_framework.predictor import Predictor
import os
import rmgpy
import logging
from rmgpy.rmg.main import initializeLog

level = logging.INFO
initializeLog(level, os.path.join('train_cnn_results', 'train0.log'))

h298_predictor = Predictor()

predictor_input = os.path.join(os.path.dirname(rmgpy.__file__),
											'cnn_framework',
											'test_data', 
											'minimal_predictor', 
											'predictor_input.py'
											)

h298_predictor.load_input(predictor_input)

h298_predictor.kfcv_train(folds=5)