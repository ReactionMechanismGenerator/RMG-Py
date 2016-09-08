#!/bin/bash
#set a job name
SBATCH --job-name=OOHAbstraction
#################
#a file for job output, you can check job progress
SBATCH --output=OOHAbstraction.out
#################
#a file for errors from the job
SBATCH --error=OOHAbstraction.err
#################
#time you think you need; default is one day
#in minutes in this case, hh:mm:ss
#SBATCH --time=20:00
#################
#number of tasts you are requesting
SBATCH -n 32
SBATCH --exclusive
#################
#partitions to use
SBATCH --partition=west
#################
# number of nodes to distribute n tasks across
SBATCH -N 2
#################

work = /gss_gpfs_scratch/harms.n/OOH_Abstraction

cd $work

mpirun -prot -srun .mpi_mm
