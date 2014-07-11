CHEMKIN MECHANISM IMPORTER INSTRUCTIONS
=======================================

To run locally, run import.sh then point a web browser at http://localhost:8084
The import.sh script assumes you have the environment variable $RMGpy pointing to 
the folder containing importChemkin.py and the rest of RMG.

If you want to run this on a cluster, then you will have to open an SSH tunnel with
port forwarding to the compute node to access the web interface. Here is an example
with the platform LSF queue system:

$ qsub import.sh
Job <18533> is submitted to default queue <medium_priority>.
Request <18533> is submitted to default queue <medium_priority>.

Then use `bjobs` to find out which compute node it's running on:

$ bjobs -w
JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
18533   myuser   RUN   medium_pri venture     compute036  import.sh  Aug  9 16:00

in this case, my job is running on compute036. Remember that compute036 is only 
accessible from the head node, not from outside the cluster. Log in to the cluster 
again, but this time with an SSH tunnel so that port 9999 on your local host is 
forwarded to port 8084 on compute036. (be sure to run this from your local computer, 
not from the cluster)

$ ssh myuser@my.cluster.edu -L 9999:compute036:8084

Then point a web browser at http://localhost:9999


To get started importing a new model, you will need to create a SMILES.txt file
and put in a few small molecules by hand:

 1. A few special-cases of things that often have different electronic states, but 
    that have the same SMILES string. For example, you might put:
    
    H2CC    singletC=[C]
    CH2     triplet[CH2]
    CH2*    singlet[CH2]
    
    These are hard-coded exceptions and can be found in the importChemkin.py source code.

 2. Any molecules which have been used as third-body colliders in pressure-dependent
    reactions in the chemkin file, but which are not strictly unambiguous from the 
    chemical formula. The most common example is:
    
    C2H2     C#C

You will most probably also have to fix a bunch of errors in the input and thermo 
files, because it seems we are more strict in reading these than chemkin is, and 
most published models contain errors.

You may also have to comment out a few reactions that RMG is unable to read:

 1. Anything with an explicit third body collider, such as:

    H+O2(+AR) <=> OOH(+AR) 1.475e+12 0.600 0.000e+00
    LOW/ 6.810e+18 -1.200 0.0/
    TROE/ 0.70 1.0e-30 1.0e+30 1.0e+30/

 2. Anything that emits a photon, such as:
 
    OH* <=> OH + Hv