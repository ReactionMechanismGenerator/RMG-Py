#!/bin/bash
#####################i################################################
# This is a job submission file for a SGE queuing system to run
# the SCOOP-enabled parallel version of RMG-Py across 48 CPUs on
# a single node.
#
# Define RMGPy as the path to rmg.py in your ~/.bash_profile
# NSLOTS is an SGE env. variable for total number of CPUs.
# prolog.sh is a script used by SCOOP to pass env. variables
#
# You can run the jobs on different nodes as well, but it is not
# recommended since you might have problems with SGE job termination.
# Type `qconf -spl` to see available parallel environments and modify
# the last #$ line if you really want to run it on many nodes.
#####################i################################################
#$ -S /bin/bash
#$ -cwd
#$ -notify
#$ -o job.log -j y
#$ -N RMGscoop
#$ -l normal
#$ -l h_rt=09:05:00 
#$ -pe singlenode 48
source ~/.bash_profile
python -m scoop.__main__ --tunnel --prolog $RMGpy/examples/rmg/scoop/prolog.sh  -n $NSLOTS $RMGpy/rmg.py input.py > std.out
