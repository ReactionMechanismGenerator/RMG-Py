# install_rms.sh
#
# Convenience script to install ReactionMechanismSimulator inside your rmg_env conda environment
#
# Note that you will have to manually install juliaup before running this script:
# curl -fsSL https://install.julialang.org | sh
# # restart shell
# juliaup add 1.10
# juliaup default 1.10
# juliaup remove release

#!/bin/bash

conda install -y conda-forge::pyjuliacall

export JULIA_CONDAPKG_BACKEND=Null
export JULIA_PYTHONCALL_EXE=$CONDA_PREFIX/bin/python

julia -e 'using Pkg; Pkg.add(Pkg.PackageSpec(name="ReactionMechanismSimulator", url="https://github.com/hwpang/ReactionMechanismSimulator.jl.git", rev="fix_installation")); using ReactionMechanismSimulator; Pkg.instantiate()' || echo "RMS install error - continuing anyway ¯\_(ツ)_/¯"

export PYTHON_JULIAPKG_EXE=$(which julia)
export PYTHON_JULIAPKG_PROJECT=$HOME/.julia/packages

echo "Copy the following text into your terminal profile (.bashrc, .zshrc, etc.):

export JULIA_CONDAPKG_BACKEND=Null
export JULIA_PYTHONCALL_EXE=$CONDA_PREFIX/bin/python
export PYTHON_JULIAPKG_EXE=$(which julia)
export PYTHON_JULIAPKG_PROJECT=$HOME/.julia/packages

or otherwise run these 4 commands when first opening a terminal that will run RMG requiring RMS.
"""
