***************
Running MEASURE
***************

Calculating Rate Coefficients
=============================

The primary mode of MEASURE calculates the phenomenological rate coefficients
for a reaction network specified in an input file. To do so, run Python on the
main ``measure.py`` file, passing the input file as a parameter on the command
line, e.g. ::

    $ python measure.py FILE

When MEASURE execution is complete, it will save an output file in the same
directory as the input file but with the name ``output.py``. If you wish to
change the location of the output, use, the ``-o``/``--output`` flag, e.g. ::

    $ python measure.py INPUTFILE -o OUTPUTFILE

Drawing Potential Energy Surface
================================

MEASURE contains functionality for automatically generating an image of the
potential energy surface for a reaction network. To do so, run Python on the
main ``measure.py`` file, passing the input file as a parameter on the command
line, and also add the ``-d``/``--draw`` flag along with the output image, 
e.g. ::

    $ python measure.py INPUTFILE -d OUTPUTIMAGE

The potential energy surface can be generated in PNG, SVG, PDF, and PS formats.
The format is automatically chosen based on the output image file extension.
If your input file contains structure information (e.g. SMILES or InChI 
strings), then the generated image will use structures to depict the species;
otherwise it will use the string labels.

Other Options
=============

Log Verbosity
-------------

You can manipulate the amount of information logged to the console window using
the ``-q``/``--quiet`` flag (for quiet mode) or the ``-v``/``--verbose`` flag 
(for verbose mode). The former causes the amount of logging information shown 
to decrease; the latter causes it to increase.

Help
----

To view help information and all available options, use the ``-h``/``--help`` 
flag, e.g. ::

    $ python measure.py -h
