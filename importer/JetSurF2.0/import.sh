#!/bin/bash
#BSUB -J port8084
#BSUB -oo output.log
#BSUB -eo error.log
#BSUB -n 1

# Run under a profiler
# Serve on Port 8084
python -m cProfile -o importChemkin.profile $rmg/importChemkin.py \
	--species Mech_JetSurF2.0.txt \
	--reactions Mech_JetSurF2.0.txt \
	--thermo Thermdat.txt \
	--known JetSurF_v1.1_SMILES.txt \
	--port 8084 --noqm
gprof2dot -f pstats  importChemkin.profile | dot -Tpdf -o importChemkin.profile.pdf
