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

make install
export PYTHON=$PREFIX/bin/python
export PYTHONPATH=$SRC_DIR:$PYTHONPATH
python -c "import julia; julia.install(); import diffeqpy; diffeqpy.install()"
julia -e 'using Pkg; Pkg.add(PackageSpec(name="ReactionMechanismSimulator",rev="main")); using ReactionMechanismSimulator'

set +x
