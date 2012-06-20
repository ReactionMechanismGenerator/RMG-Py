import os

RMG_path = os.environ.get("RMGpy")
if RMG_path is None:
    RMG_path = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..'))
    import logging
    logging.info("Setting RMG_path to {0}".format(RMG_path))

RMG_bin_path = os.path.join(RMG_path, 'bin')
if not os.path.exists(RMG_bin_path):
    raise Exception("Directory {0} does not exist.".format(RMG_bin_path))
if not os.path.isdir(RMG_bin_path):
    raise Exception("{0} is not a directory.".format(RMG_bin_path))

