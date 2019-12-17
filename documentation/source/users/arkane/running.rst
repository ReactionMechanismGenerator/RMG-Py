**************
Running Arkane
**************

To execute an Arkane job, invoke the command ::

    $ python Arkane.py INPUTFILE

The absolute or relative paths to the Arkane.py file as well as to the input file must be given.

The job will run and the results will be saved in the same
directory as the input file. If you wish to save the output elsewhere, use
the ``-o``/``--output`` option, e.g. ::

    $ python Arkane.py INPUTFILE -o OUTPUTFILE

Drawing Potential Energy Surface
================================

Arkane contains functionality for automatically generating an image of the
potential energy surface for a reaction network. This is done automatically
and outputted in pdf format to a file called ``network.pdf``.


Log Verbosity
=============

You can manipulate the amount of information logged to the console window using
the ``-q``/``--quiet`` flag (for quiet mode), the ``-v``/``--verbose`` flag
(for verbose mode), or the ``-d``/``--debug`` flag (for debug mode).
The former causes the amount of logging information shown
to decrease; the latter two cause it to increase.

Suppressing plot creation
=========================

Arkane by default will generate many plot files. By adding the ``-p``/``--no-plot``
flag, Arkane will not generate any plots, reducing file size of output and
increasing the calculation speed.

Help
====

To view help information and all available options, use the ``-h``/``--help`` 
flag, e.g. ::

    $ python Arkane.py -h

