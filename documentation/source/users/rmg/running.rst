.. _running:

*************
Running a Job
*************

Running a basic RMG job is straightforward. However, depending on your case you might want to add the flags outlined in the following examples.

**Note:** In all these examples ``rmg.py`` should be the path to your installed RMG (eg. yours might be ``/Users/joeblogs/Code/RMG-Py/rmg.py``) and ``input.py`` is the path to the input file you wish to run (eg. yours might be ``RMG-runs/hexadiene/input.py``).  If you get an error like ``python: can't open file 'rmg.py': [Errno 2] No such file or directory``  then probably the first of these is wrong. If you get an error like ``IOError: [Errno 2] No such file or directory: '/some/path/to/input.py'`` then probably the second of these is wrong.

Basic run::

	python rmg.py input.py

Run by restarting from a seed mechanism ::

    python rmg.py input.py -r path/to/seed/

Run with CPU profiling::

    python rmg.py input.py -p

We recommend you make a job-specific directory for each RMG simulation. Some jobs can take quite a while to complete, so we also recommend using a job scheduler if working in an linux environment. 

The instructions below describe special cases for running an RMG job.

Running RMG in parallel with SLURM
----------------------------------

RMG has the option to use multiple processes on one node for reaction generation and on-the-fly Quantum Mechanics Thermodynamic Property (QMTP) calculation. Here is an example submission script for an RMG-Py job with a SLURM scheduler.

The job reserves 24 tasks on a single node, but uses only 12 processes in parallel during
the RMG-Py simulation.

Make sure that: 

- the queue named ``debug`` exists on your SLURM scheduler. 
- you modify the path to the parent folder of the RMG-Py installation folder.
- you have an anaconda environment named ``rmg_env`` that contains RMG-Py's dependencies.
- the working directory from which you launched the job contains the RMG-Py input file ``input.py``

.. code:: bash

    #!/bin/bash

    #SBATCH -p debug
    #SBATCH -J jobname
    #SBATCH -n 24

    Processes=12
    RMG_WS=/path/to/RMG/parent/folder
    export PYTHONPATH=$PYTHONPATH:$RMG_WS/RMG-Py/

    source activate rmg_env

    python $RMG_WS/RMG-Py/rmg.py -n $Processes input.py

    source deactivate

Running RMG in parallel with SGE
--------------------------------

RMG has the option to use multiple processes on one node for reaction generation and on-the-fly Quantum Mechanics Thermodynamic Property (QMTP) calculation. Here is an example submission script for an RMG-Py job with a SGE scheduler.

The job reserves 24 tasks on a single node, but uses only 12 processes in parallel during
the RMG-Py simulation.

Make sure that:

-  the queue named ``debug`` exists on your SGE scheduler.
-  you modify the path to the parent folder of the RMG-Py installation
   folder.
-  you have an anaconda environment named ``rmg_env`` that contains
   RMG-Py's dependencies.
-  the working directory from which you launched the job
   contains the RMG-Py input file ``input.py``.

.. code:: bash

    #! /bin/bash

    #$ -l debug
    #$ -N jobname
    #$ -pe singlenode 24

    Processes=12
    RMG_WS=/path/to/RMG/parent/folder
    export PYTHONPATH=$PYTHONPATH:$RMG_WS/RMG-Py/

    source activate rmg_env

    python $RMG_WS/RMG-Py/rmg.py -n $Processes input.py

    source deactivate


Details on the implementation
--------------------------------

Currently, multiprocessing is implemented for reaction generation and the generation of QMfiles when using the QMTP option to compute thermodynamic properties of species. The processes are spawned and closed within each function. The number of processes is determined based on the ratio of currently available RAM and currently used RAM. The user can input the maximum number of allowed processes from the command line. For each reaction generation or QMTP call the number of processes will be the minimum value of either the number of allowed processes due to user input or the value obtained by the RAM ratio. The RAM limitation is employed, because multiprocessing is forking the base process and the memory limit (SWAP + RAM) might be exceeded when using too many processors for a base process large in memory.

In python 3.4 new forking contexts 'spawn' and 'forkserver' are available. These methods will create new processes which share nothing or limited state with the parent and all memory passing is explicit. Once RMG is transferred to python 3 it is recommended to use the spawn or forkserver forking context to potentially allow for an increased number of processes.
