#!/bin/bash
#The script that each workers of scoop runs.
#Should be modified based on user profile.
source ~/.bashrc
export RMGQM="/opt/rmgqm"
#set Gaussian03 environment variables
g03root=/opt
GAUSS_SCRDIR=/scratch/$USER
export g03root GAUSS_SCRDIR
GAUSS_EXEDIR="$g03root/g03/"
GAUSS_LEXEDIR="$g03root/g03/linda-exe"
GAUSS_ARCHDIR="$g03root/g03/arch"
GMAIN=$GAUSS_EXEDIR/g03
PATH=$PATH:$GMAIN
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GMAIN
G03BASIS="$g03root/g03/basis"
F_ERROPT1="271,271,2,1,2,2,2,2"
TRAP_FPE="OVERFL=ABORT;DIVZERO=ABORT;INT_OVERFL=ABORT"
MP_STACK_OVERFLOW="OFF"
# to partially avoid KAI stupidity
KMP_DUPLICATE_LIB_OK="TRUE"
export GAUSS_EXEDIR GAUSS_ARCHDIR PATH GMAIN LD_LIBRARY_PATH F_ERROPT1 TRAP_FPE MP_STACK_OVERFLOW \
  KMP_DUPLICATE_LIB_OK G03BASIS GAUSS_LEXEDIR
#set MOPAC
export MOPAC_LICENSE=/opt/mopac/
#
export PYTHONPATH=/home/keceli/RMG/RMG-Py/PyDAS/build/lib.linux-x86_64-2.6:/home/keceli/RMG/RMG-Py/PyDQED:/home/keceli/local/lib/python2.6/site-packages:/opt/rmgqm/RDKit_2013_03_2:/usr/local/lib/python2.6/dist-packages:$PYTHONPATH
export PATH=/home/keceli/kiler:/home/keceli/bin:/home/keceli/local/bin:/opt/mpich2-1.2.1p1/bin:/opt/intel/Compiler/11.0/074/bin/intel64:/opt/sge/bin/lx24-amd64:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/bin/mh:/opt/g03/g03:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/Compiler/11.0/074/ipp/em64t/sharedlib:/opt/intel/Compiler/11.0/074/mkl/lib/em64t:/opt/intel/Compiler/11.0/074/tbb/em64t/cc4.1.0_libc2.4_kernel2.6.16.21/lib:/opt/intel/Compiler/11.0/074/lib/intel64:/opt/g03/g03:/usr/local/lib:/opt/rmgqm/RDKit_2013_03_2/bin:/opt/rmgqm/boost_1_44_0/lib:/opt/rmgqm/RDKit_2013_03_2/lib
