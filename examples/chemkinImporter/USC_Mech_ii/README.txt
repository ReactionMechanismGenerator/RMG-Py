To run locally, run import.sh then point a web browser at http://localhost:8084

If you want to run this on a cluster, then you will have to open an SSH tunnel with port forwarding to the compute node to access the web interface. Here is an example with the platform LSF queue system:

$ qsub import.sh
Job <18533> is submitted to default queue <medium_priority>.
Request <18533> is submitted to default queue <medium_priority>.

Then use `bjobs` to find out which compute node it's running on:

$ bjobs
JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
18533   myuser   RUN   medium_pri venture     compute036  import.sh  Aug  9 16:00

in this case, my job is running on compute036. Remember that compute036 is only accessible from the head node, not from outside the cluster. Log in to the cluster again, but this time with an SSH tunnel so that port 9999 on your local host is forwarded to port 8084 on compute036. (be sure to run this from your local computer, not from the cluster)

$ ssh myuser@my.cluster.edu -L 9999:compute036:8084

Then point a web browser at http://localhost:9999

