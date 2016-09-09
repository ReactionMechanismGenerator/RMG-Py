#!/bin/sh
#set a job name
#SBATCH --job-name=AutoTST

#a file for job output, you can check job progress
#SBATCH --output=AutoTST.out.%a.log

#a file for errors from the job
#SBATCH --error=AutoTST.error.%a.log

#time you think you need; default is one day
#hh:mm:ss
#SBATCH --time=05:00:00

#number of tasks you are requesting
#SBATCH -n 11
#SBATCH --exclusive

#partition to use
#SBATCH --partition=west

#number of nodes to distribute n tasks across
#SBATCH -N 1

#an array job
#SBATCH --array=1-20


#####################################################

export RMGpy=~/Code/RMG-Py
export PYTHONPATH=$RMGpy:$PYTHONPATH
cd $RMGpy/scripts
#python filterReactions.py /scratch/westgroup/Importer/RMG-models/
## that creates the kineticsDict files, and doesn't need repeating until the imported models change significantly
echo $SLURM_ARRAY_TASK_ID
python autoTST-OOH.py


