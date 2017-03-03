
import os
import logging
from .cnn_model import build_model

predictor = None

def predictor_model(embedding_size=512, attribute_vector_size=None, depth=2, 
                scale_output=0.05, padding=False, 
                hidden=0, hidden_activation='tanh',
                output_activation='linear', output_size=1, 
                lr=0.01, optimizer='adam', loss='mse'):
    
    model = build_model(embedding_size, attribute_vector_size, depth, 
                scale_output, padding, 
                hidden, hidden_activation,
                output_activation, output_size, 
                lr, optimizer, loss)
    
    predictor.model = model

def read_input_file(path, predictor0):

    global predictor
    
    full_path = os.path.abspath(os.path.expandvars(path))
    try:
        f = open(full_path)
    except IOError, e:
        logging.error('The input file "{0}" could not be opened.'.format(full_path))
        logging.info('Check that the file exists and that you have read access.')
        raise e

    logging.info('Reading predictor input file "{0}"...'.format(full_path))
    logging.info(f.read())
    f.seek(0)# return to beginning of file

    predictor = predictor0
    
    global_context = { '__builtins__': None }
    local_context = {
        '__builtins__': None,
        'True': True,
        'False': False,
        'predictor_model': predictor_model,
    }

    try:
        exec f in global_context, local_context
    except (NameError, TypeError, SyntaxError), e:
        logging.error('The input file "{0}" was invalid:'.format(full_path))
        logging.exception(e)
        raise
    finally:
        f.close()

    logging.info('')