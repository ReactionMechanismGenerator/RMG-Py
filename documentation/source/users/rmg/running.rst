.. _running:

*************
Running a Job
*************

**Note:** In all these examples ``rmg.py`` should be the path to your installed RMG (eg. yours might be ``/Users/joeblogs/Code/RMG-Py/rmg.py``) and ``input.py`` is the path to the input file you wish to run (eg. yours might be ``RMG-runs/hexadiene/input.py``).  If you get an error like ``python: can't open file 'rmg.py': [Errno 2] No such file or directory``  then probably the first of these is wrong. If you get an error like ``IOError: [Errno 2] No such file or directory: '/some/path/to/input.py'`` then probably the second of these is wrong.

Running a basic RMG job is straightforward, as shown in the example below. However, depending on your case you might want to add the flags outlined in the following section. We recommend you make a job-specific directory for each RMG simulation. Some jobs can take quite a while to complete, so we also recommend using a job scheduler if working in a linux environment. 

Basic run::

	python rmg.py input.py

.. _inputflags:

Input flags
-----------

The options for input flags can be found in ``/RMG-Py/rmgpy/util.py``. Running ::

 	python rmg.py -h

at the command line will print the documentation from ``util.py``, which is reproduced below for convenience::

	usage: rmg.py [-h] [-q | -v | -d] [-o DIR] [-r path/to/seed/] [-p] [-P]
              [-t DD:HH:MM:SS] [-i MAXITER] [-n MAXPROC] [-k]
              FILE

	Reaction Mechanism Generator (RMG) is an automatic chemical reaction mechanism
	generator that constructs kinetic models composed of elementary chemical
	reaction steps using a general understanding of how molecules react.

	positional arguments:
	  FILE                  a file describing the job to execute

	optional arguments:
	  -h, --help            show this help message and exit
	  -q, --quiet           only print warnings and errors
	  -v, --verbose         print more verbose output
	  -d, --debug           print debug information
	  -o DIR, --output-directory DIR
	                        use DIR as output directory
	  -r path/to/seed/, --restart path/to/seed/
	                        restart RMG from a seed
	  -p, --profile         run under cProfile to gather profiling statistics, and
	                        postprocess them if job completes
	  -P, --postprocess     postprocess profiling statistics from previous
	                        [failed] run; does not run the simulation
	  -t DD:HH:MM:SS, --walltime DD:HH:MM:SS
	                        set the maximum execution time
	  -i MAXITER, --maxiter MAXITER
	                        set the maximum number of RMG iterations
	  -n MAXPROC, --maxproc MAXPROC
	                        max number of processes used during reaction
	                        generation
	  -k, --kineticsdatastore
	                        output a folder, kinetics_database, that contains a
	                        .txt file for each reaction family listing the
	                        source(s) for each entry

Some representative example usages are shown below.

Run by restarting from a seed mechanism::

    python rmg.py -r path/to/seed/ input.py

Run with CPU time profiling::

    python rmg.py -p input.py

Run with multiprocessing for reaction generation and QMTP::

    python rmg.py -n <Max number of processes allowed> input.py 

Run with setting a limit on the maximum execution time::

	python rmg.py -t <DD:HH:MM:SS> input.py

Run with setting a limit on the maximum number of iterations::

	python rmg.py -i <Max number of desired iterations> input.py


Details on the multiprocessing implementation
---------------------------------------------

Currently, multiprocessing is implemented for reaction generation and the generation of QMfiles when using the QMTP option to compute thermodynamic properties of species. The processes are spawned and closed within each function. The number of processes is determined based on the ratio of currently available RAM and currently used RAM. The user can input the maximum number of allowed processes from the command line. For each reaction generation or QMTP call the number of processes will be the minimum value of either the number of allowed processes due to user input or the value obtained by the RAM ratio. The RAM limitation is employed, because multiprocessing is forking the base process and the memory limit (SWAP + RAM) might be exceeded when using too many processors for a base process large in memory.

In python 3.4 new forking contexts 'spawn' and 'forkserver' are available. These methods will create new processes which share nothing or limited state with the parent and all memory passing is explicit. Once RMG is transferred to python 3 it is recommended to use the spawn or forkserver forking context to potentially allow for an increased number of processes.

