set -x

ACTIVATE_ENV="${PREFIX}/etc/conda/activate.d/env_vars.sh"

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

# make install -> should we 'install' or 'all'?
make
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
