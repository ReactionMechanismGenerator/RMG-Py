#!/bin/sh
#BSUB -o RMG.out
#BSUB -J RMGPyScoop
#BSUB -n 8
#BSUB -e error_log
#BSUB -q medium_priority

# This is a job submission file for a LSF queuing system to run
# the SCOOP-enabled parallel version of RMG-Py across 8 CPUs on
# a number of different compute nodes on a (potentially heterogeneous) cluster.

source ~/.bash_profile

LAMHOST_FILE=hosts

# start a new host file from scratch
rm -f $LAMHOST_FILE
touch $LAMHOST_FILE
# echo "# LAMMPI host file created by LSF on `date`" >> $LAMHOST_FILE
# check if we were able to start writing the conf file
if [ -f $LAMHOST_FILE ]; then
	:
else
	echo "$0: can't create $LAMHOST_FILE"
	exit 1
fi
HOST=""
NUM_PROC=""
FLAG=""
TOTAL_CPUS=0
for TOKEN in $LSB_MCPU_HOSTS
do
	if [ -z "$FLAG" ]; then
		HOST="$TOKEN"
		FLAG="0"
	else
		NUM_PROC="$TOKEN"
		TOTAL_CPUS=`expr $TOTAL_CPUS + $NUM_PROC`
		FLAG="1"
	fi
	if [ "$FLAG" = "1" ]; then
		_x=0
		while [ $_x -lt $NUM_PROC ]
		do
			echo "$HOST" >>$LAMHOST_FILE
			_x=`expr $_x + 1`
		done
		# get ready for the next host
		FLAG=""
		HOST=""
		NUM_PROC=""
	fi
done
# last thing added to LAMHOST_FILE
#echo "# end of LAMHOST file" >> $LAMHOST_FILE
echo "Your lamboot hostfile looks like:"
cat $LAMHOST_FILE

python -m scoop -vv --hostfile $LAMHOST_FILE $RMGpy/rmg.py input.py > RMG.stdout.log
