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

### Julia breaks the Cantera installation (by changing the installed version of Sundials) - this prevents that (hopefully)
# Detect platform
OS=$(uname)
IS_MAC=false
if [ "$OS" = "Darwin" ]; then
    IS_MAC=true
fi
# Determine extension based on platform
if [ "$IS_MAC" = true ]; then
    LIB_EXT="dylib"
else
    LIB_EXT="so"
fi
# Find the Cantera-linked shared library
SUNDIALS_SO_PATH=$(find "$CONDA_PREFIX" -name "libcantera*.${LIB_EXT}*" 2>/dev/null | head -n 1)
if [ -z "$SUNDIALS_SO_PATH" ]; then
    echo "❌ Cantera-linked SUNDIALS library not found in the current conda environment."
    exit 1
fi
SUNDIALS_SO_DIR=$(dirname "$SUNDIALS_SO_PATH")
# Export the correct environment variable
# the sundials docs specifically request setting this environment variable
# (https://github.com/SciML/Sundials.jl?tab=readme-ov-file#installation)
# rather than the more common LD_LIBRARY_PATH
# for macos we set DYLD_LIBRARY_PATH, which is (allegedly) the Mac equivalent
# though this might be the wrong thing to do
if [ "$IS_MAC" = true ]; then
    export DYLD_LIBRARY_PATH="$SUNDIALS_SO_DIR:$DYLD_LIBRARY_PATH"
    echo "✅ DYLD_LIBRARY_PATH updated for macOS."
else
    export DL_LOAD_PATH="$SUNDIALS_SO_DIR:$DL_LOAD_PATH"
    echo "✅ DL_LOAD_PATH updated for Linux."
fi
###

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
