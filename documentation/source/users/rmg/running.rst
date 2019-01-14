.. _running:

*************
Running a Job
*************

Running RMG job is easy and under different situations you might want add additional flag as the following examples.

**Note:** In all these examples ``rmg.py`` should be the path to your installed RMG (eg. yours might be ``/Users/joeblogs/Code/RMG-Py/rmg.py``) and ``input.py`` is the path to the input file you wish to run (eg. yours might be ``RMG-runs/hexadiene/input.py``).  If you get an error like ``python: can't open file 'rmg.py': [Errno 2] No such file or directory``  then probably the first of these is wrong. If you get an error like ``IOError: [Errno 2] No such file or directory: '/some/path/to/input.py'`` then probably the second of these is wrong.

Basic run::

	python rmg.py input.py

Run with a restart file (restart file should be located in same folder as input.py)::

    python rmg.py input.py -r

Run with CPU profiling::

    python rmg.py input.py -p

We recommend you make a job-specific directory for each RMG simulation. Some jobs can take quite a while to complete, so we also recommend using a job scheduler (if working in an linux environment). 

The instructions below describe more special cases for running an RMG job.

Running RMG in parallel with SLURM
----------------------------------

RMG has the capability to run using multiple cores. Here is an example
job submission script for an RMG-Py job with a SLURM scheduler

The job named ``min_par`` reserves 24 CPUs on a single node
(``-np 24``), but uses only 12 workers (= 12 CPUs) in parallel during
the RMG-Py simulation.

Make sure that: 

- the queue named ``debug`` exists on your SLURM scheduler. 
- you modify the path to the parent folder of the RMG-Py installation folder 
- you have an anaconda environment named ``rmg_env`` that contains RMG-Py's dependencies 
- the working directory from which you launched the job contains the RMG-Py input file ``input.py``


``-v`` adds verbosity to the output log file.

.. code:: bash

    #!/bin/bash
    #SBATCH -p debug
    #SBATCH -J min_par
    #SBATCH -n 24

    hosts=$(srun bash -c hostname)

    WORKERS=12

    RMG_WS=/path/to/RMG/parent/folder
    export PYTHONPATH=$PYTHONPATH:$RMG_WS/RMG-Py/

    source activate rmg_env
    python -m scoop -n $WORKERS --host $hosts -v $RMG_WS/RMG-Py/rmg.py input.py
    source deactivate

Running RMG in parallel with SGE
--------------------------------

RMG has the capability to run using multiple cores. Here is an example
using the SGE scheduler.

In order to help understand, the example job is also named ``min_par``
reserving 24 CPUs on a single node (``#$ -pe singlenode 24``), but uses
only 12 workers (= 12 CPUs) in parallel during the RMG-Py simulation.

Make sure that:

-  the queue named ``normal`` exists on your SGE scheduler
-  you modify the path to the parent folder of the RMG-Py installation
   folder
-  you have an anaconda environment named ``rmg_env`` that contains
   RMG-Py's dependencies
-  the working directory from which you launched the job
   contains the RMG-Py input file ``input.py``

``-v`` adds verbosity to the output log file

.. code:: bash

    #! /bin/bash

    #$ -o job.log
    #$ -l normal
    #$ -N min_par
    #$ -pe singlenode 24

    WORKERS=12

    RMG_WS=/path/to/RMG/parent/folder
    export PYTHONPATH=$PYTHONPATH:$RMG_WS/RMG-Py/

    source activate rmg_env
    python -m scoop --tunnel -n $WORKERS -v $RMG_WS/RMG-Py/rmg.py input.py

    source deactivate

