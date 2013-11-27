#!/bin/bash
# Run under a profiler
# Serve on Port 8084

python -m cProfile -o importChemkin.profile $RMGpy/importChemkin.py  --species USC_Mech_ver_II.txt --reactions USC_Mech_ver_II.txt --thermo thermdat.txt --known SMILES.txt --port 8084
gprof2dot -f pstats  importChemkin.profile | dot -Tpdf -o importChemkin.profile.pdf