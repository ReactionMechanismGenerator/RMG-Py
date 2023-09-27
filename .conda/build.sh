# conda build script
#
# conda executes this script during the build process to make the rmg
# binary
make install

# The below code does not work because the Julia programming langauge is a very interesting technical demo that is completely unfit for actual
# deployment in any real world projects.
# 
# RMG will be shipped as a pure python package, and then RMS installed by the user.
#
# export PYTHON=$PREFIX/bin/python
# export PYTHONPATH=$SRC_DIR:$PYTHONPATH
# python -c "import julia; julia.install(); import diffeqpy; diffeqpy.install()"
# julia -e 'using Pkg; Pkg.add(PackageSpec(name="Functors",version="0.4.3")); Pkg.pin("Functors"); Pkg.add(PackageSpec(name="ReactionMechanismSimulator",rev="for_rmg")); using ReactionMechanismSimulator'
