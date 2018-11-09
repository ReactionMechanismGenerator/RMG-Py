**************
Running Arkane
**************

To execute an Arkane job, invoke the command ::

    $ python arkane.py INPUTFILE

The absolute or relative paths to the arkane.py file as well as to the input file must be given.

The job will run and the results will be saved to ``output.py`` in the same
directory as the input file. If you wish to save the output elsewhere, use
the ``-o``/``--output`` option, e.g. ::

    $ python arkane.py INPUTFILE -o OUTPUTFILE

Drawing Potential Energy Surface
================================

Arkane contains functionality for automatically generating an image of the
potential energy surface for a reaction network. This is done automatically
and outputted in pdf format to a file called ``network.pdf``.


Log Verbosity
=============

You can manipulate the amount of information logged to the console window using
the ``-q``/``--quiet`` flag (for quiet mode) or the ``-v``/``--verbose`` flag
(for verbose mode). The former causes the amount of logging information shown
to decrease; the latter causes it to increase.

Help
====

To view help information and all available options, use the ``-h``/``--help`` 
flag, e.g. ::

    $ python arkane.py -h

