
from rmgpy.cnn_framework.predictor import Predictor
import os
import rmgpy
import logging
import argparse

def parseCommandLineArguments():
    """
    Parse the command-line arguments being passed to RMG Py. This uses the
    :mod:`argparse` module, which ensures that the command-line arguments are
    sensible, parses them, and returns them.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('file', metavar='FILE', type=str, nargs=1,
        help='a predictor training input file')

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

	input_file = args.file[0]
	input_directory = os.path.abspath(os.path.dirname(input_file))
	
	level = logging.INFO
	initializeLog(level, os.path.join(input_directory, 'train.log'))

	h298_predictor = Predictor()

	h298_predictor.load_input(input_file)

	lr_func = "float(0.0007 * np.exp(- epoch / 30.0))"
	h298_predictor.kfcv_train(folds=5, lr_func=lr_func)