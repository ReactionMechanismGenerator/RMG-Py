
from rmgpy.cnn_framework.predictor import Predictor
import os
import sys
import rmgpy
import logging
import argparse
import shutil
import time

def parseCommandLineArguments():
	"""
	Parse the command-line arguments being passed to RMG Py. This uses the
	:mod:`argparse` module, which ensures that the command-line arguments are
	sensible, parses them, and returns them.
	"""

	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', metavar='FILE', type=str, 
		nargs=1, help='a predictor training input file')

	parser.add_argument('-d', '--datasets', metavar='FILE', type=str, 
		nargs='+', help='a file specifies on which datasets to train')

	parser.add_argument('-f', '--folds', type=int, 
		default=5, help='number of folds for training')

	parser.add_argument('-t', '--train_mode', type=str, 
		help='train mode: currently support in_house and keras')

	parser.add_argument('-bs', '--batch_size', type=int, 
		help='batch training size')

	parser.add_argument('-lr', '--learning_rate', type=str, 
		default='0.0007_30.0', help='two parameters for learning rate')

	return parser.parse_args()
################################################################################

def initializeLog(verbose, log_file_name):
	"""
	Set up a logger to print output to stdout. The
	`verbose` parameter is an integer specifying the amount of log text seen
	at the console; the levels correspond to those of the :data:`logging` module.
	"""
	# Create logger
	logger = logging.getLogger()
	logger.setLevel(verbose)

	# Create console handler and set level to debug; send everything to stdout
	# rather than stderr
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(verbose)

	logging.addLevelName(logging.CRITICAL, 'Critical: ')
	logging.addLevelName(logging.ERROR, 'Error: ')
	logging.addLevelName(logging.WARNING, 'Warning: ')
	logging.addLevelName(logging.INFO, '')
	logging.addLevelName(logging.DEBUG, '')
	logging.addLevelName(1, '')

	# Create formatter and add to console handler
	#formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', '%Y-%m-%d %H:%M:%S')
	#formatter = Formatter('%(message)s', '%Y-%m-%d %H:%M:%S')
	formatter = logging.Formatter('%(levelname)s%(message)s')
	ch.setFormatter(formatter)

	# create file handler
	if os.path.exists(log_file_name):
		backup = os.path.join(log_file_name[:-9], 'train_backup.log')
		if os.path.exists(backup):
			print "Removing old "+backup
			os.remove(backup)
		print 'Moving {0} to {1}\n'.format(log_file_name, backup)
		shutil.move(log_file_name, backup)
	fh = logging.FileHandler(filename=log_file_name) #, backupCount=3)
	fh.setLevel(min(logging.DEBUG,verbose)) # always at least VERBOSE in the file
	fh.setFormatter(formatter)
	# notice that STDERR does not get saved to the log file
	# so errors from underlying libraries (eg. openbabel) etc. that report
	# on stderr will not be logged to disk.

	# remove old handlers!
	while logger.handlers:
		logger.removeHandler(logger.handlers[0])

	# Add console and file handlers to logger
	logger.addHandler(ch)
	logger.addHandler(fh)

################################################################################

if __name__ == '__main__':

	# to run the script
	# example command: 
	# python train_cnn.py -i input.py -d datasets.txt -f 5 -t in_house -bs 1 -lr 0.0007_30.0

	args = parseCommandLineArguments()
	input_file = args.input[0]
	datasets_file = args.datasets[0]
	folds = args.folds
	train_mode = args.train_mode
	batch_size = args.batch_size
	lr0, lr1 = [float(i) for i in args.learning_rate.split('_')]

	input_directory = os.path.abspath(os.path.dirname(input_file))
	
	level = logging.INFO
	initializeLog(level, os.path.join(input_directory, 'train.log'))

	# Log start timestamp
	logging.info('CNN training initiated at ' + time.asctime() + '\n')

	from rmgpy.rmg.main import RMG
	rmg = RMG()
	rmg.logHeader()

	predictor = Predictor(datasets_file=datasets_file)
	predictor.load_input(input_file)

	lr_func = "float({0} * np.exp(- epoch / {1}))".format(lr0, lr1)
	save_model_path = os.path.join(input_directory, 'saved_model')
	if not os.path.exists(save_model_path):
		os.mkdir(save_model_path)
	
	if train_mode == 'in_house':
		predictor.kfcv_train(folds=folds, 
							 batch_size=batch_size, 
							 lr_func=lr_func, 
							 save_model_path=save_model_path)
	elif train_mode == 'keras':
		predictor.kfcv_batch_train(folds=folds, batch_size=batch_size)
	elif train_mode == 'full_train':
		predictor.full_train(batch_size=batch_size, 
							 lr_func=lr_func, 
							 save_model_path=save_model_path)
	else:
		raise Exception('Currently not supporting train mode: {0}'.format(train_mode))