import os
if not os.environ.get("RMG_workingDirectory"):
    import os.path
    message = "Please set your RMG_workingDirectory environment variable.\n" +\
        "(eg. export RMG_workingDirectory=%s )" % \
        os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..'))
    raise Exception(message)
