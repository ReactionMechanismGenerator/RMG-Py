# from https://github.com/ReactionMechanismGenerator/RMG-Py/pull/2631#issuecomment-1998723914
#make julia directory
mkdir -p ${PREFIX}/share/julia/site
mkdir -p ${PREFIX}/bin
#set JULIA_DEPOT_PATH in conda env
export JULIA_DEPOT_PATH="${PREFIX}/share/julia/site"
ACTIVATE_ENV="${PREFIX}/etc/conda/activate.d/env_vars.sh"
DEACTIVATE_ENV="${PREFIX}/etc/conda/deactivate.d/env_vars.sh"

if [ -f "$ACTIVATE_ENV" ]; then
        echo "export JULIA_DEPOT_PATH=\"${PREFIX}/share/julia/site\"" >> $ACTIVATE_ENV
        echo "export JULIA_OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH" >> $ACTIVATE_ENV
        echo "export LD_LIBRARY_PATH=\"${PREFIX}/lib\"" >> $ACTIVATE_ENV
else
        mkdir -p ${PREFIX}/etc/conda/activate.d
        touch ${PREFIX}/etc/conda/activate.d/env_vars.sh
        echo '#!/bin/sh' >> $ACTIVATE_ENV
        echo "export JULIA_DEPOT_PATH=\"${PREFIX}/share/julia/site\"" >> $ACTIVATE_ENV
        echo "export JULIA_OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH" >> $ACTIVATE_ENV
        echo "export LD_LIBRARY_PATH=\"${PREFIX}/lib\"" >> $ACTIVATE_ENV
fi
if [ -f "$DEACTIVATE_ENV" ]; then
        echo "unset JULIA_DEPOT_PATH" >> $DEACTIVATE_ENV
        echo "export LD_LIBRARY_PATH=$JULIA_OLD_LD_LIBRARY_PATH" >> $DEACTIVATE_ENV
        echo "unset JULIA_OLD_LD_LIBRARY_PATH" >> $DEACTIVATE_ENV
else
        mkdir -p ${PREFIX}/etc/conda/deactivate.d
        touch ${PREFIX}/etc/conda/deactivate.d/env_vars.sh
        echo '#!/bin/sh' >> $DEACTIVATE_ENV
        echo "unset JULIA_DEPOT_PATH" >> $DEACTIVATE_ENV
        echo "export LD_LIBRARY_PATH=$JULIA_OLD_LD_LIBRARY_PATH" >> $DEACTIVATE_ENV
        echo "unset JULIA_OLD_LD_LIBRARY_PATH" >> $DEACTIVATE_ENV
fi

make install
export PYTHON=$PREFIX/bin/python
export PYTHONPATH=$SRC_DIR:$PYTHONPATH
python -c "import julia; julia.install(); import diffeqpy; diffeqpy.install()"
julia -e 'using Pkg; Pkg.add(PackageSpec(name="ReactionMechanismSimulator",rev="main")); using ReactionMechanismSimulator'