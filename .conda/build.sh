set -x
# as done in the pure-python conda recipe
make install
python -c "from rmgpy.molecule import Molecule"

# attempting to run the Julia install
# from https://github.com/ReactionMechanismGenerator/RMG-Py/pull/2631#issuecomment-1998723914
#make julia directory
mkdir -p ${PREFIX}/share/julia/site
mkdir -p ${PREFIX}/bin
#set JULIA_DEPOT_PATH in conda env
export JULIA_DEPOT_PATH="${PREFIX}/share/julia/site"
ACTIVATE_ENV="${PREFIX}/etc/conda/activate.d/env_vars.sh"
DEACTIVATE_ENV="${PREFIX}/etc/conda/deactivate.d/env_vars.sh"

if [ -f "$ACTIVATE_ENV" ]; then
        echo 'python-jl -c "import julia; julia.install()"' >> $ACTIVATE_ENV
        # echo 'sed -i \'/julia.install/d\' $ACTIVATE_ENV' >> $ACTIVATE_ENV
else
        mkdir -p ${PREFIX}/etc/conda/activate.d
        touch ${PREFIX}/etc/conda/activate.d/env_vars.sh
        echo '#!/bin/sh' >> $ACTIVATE_ENV
        echo 'python-jl -c "import julia; julia.install()"' >> $ACTIVATE_ENV
        # echo 'sed -i \'/julia.install/d\' $ACTIVATE_ENV' >> $ACTIVATE_ENV
fi

export PYTHON=$PREFIX/bin/python
export PYTHONPATH=$SRC_DIR:$PYTHONPATH
echo "testing rmgpy"
python -c "from rmgpy.molecule import Molecule"

echo "pythonpath"
echo $PYTHONPATH 
echo "python"
echo $PYTHON


python -c "import julia; julia.install(); import diffeqpy; diffeqpy.install()"
echo "current python"
which python
echo "julia python"
julia -e "using PyCall; println(PyCall.PYTHONHOME)"
echo "retest rmgpy"
python -c "from rmgpy.molecule import Molecule"

julia -e 'using Pkg; Pkg.add(PackageSpec(name="ReactionMechanismSimulator",rev="main")); using ReactionMechanismSimulator'

set +x
