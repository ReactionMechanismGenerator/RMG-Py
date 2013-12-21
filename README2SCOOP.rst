******************************************************
SCOOP enabled RMG-Py
******************************************************

RMG-Py can be run in parallel (only for the thermochemical parameter 
estimation part) using SCOOP module.
More info on SCOOP: http://code.google.com/p/scoop/

Running RMG-Py in parallel:

python -m scoop.__main__  -n 8 $RMGpy/rmddg.py input.py > RMG.sdout.log &

-n 8 specifies that you will have 8 workers. 
Set it based on the available number of processors.
For job submission scripts check examples/rmg/scoop.

Installing SCOOP:

You need the development version of SCOOP (tagged with 0.7RC2).
Download link: http://scoop.googlecode.com/archive/0.7RC2.zip

