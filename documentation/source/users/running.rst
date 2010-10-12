****************
Running CanTherm
****************

To execute a CanTherm job, invoke the command ::

    $ python cantherm.py FILE

The job will run and the results will be saved to ``output.py`` in the same
directory as the input file. If you wish to save the output elsewhere, use
the ``-o``/``--output`` option, e.g. ::

    $ python cantherm.py INPUTFILE -o OUTPUTFILE

You can manipulate the amount of information logged to the console window using
the ``-q``/``--quiet`` flag (for quiet mode) or the ``-v``/``--verbose`` flag 
(for verbose mode). The former causes the amount of logging information shown 
to decrease; the latter causes it to increase.

To view help information and all available options, use the ``-h``/``--help`` 
flag, e.g. ::

    $ python cantherm.py -h
