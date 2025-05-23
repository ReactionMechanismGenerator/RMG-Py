# install_rms.sh
#
# Convenience script to install ReactionMechanismSimulator into an rmg_env conda environment
#

# Check if juliaup is installed
if ! command -v juliaup &> /dev/null; then
    echo "Could not find julia via juliaup. Please install it by running:"
    echo "curl -fsSL https://install.julialang.org | sh"
    echo "# Restart your shell"
    echo "juliaup add 1.10"
    echo "juliaup default 1.10"
    echo "juliaup remove release"
    exit 1
fi

# Check if julia command is available
if ! command -v julia &> /dev/null; then
    echo "Julia is not installed. Please install version 1.10 by running:"
    echo "juliaup add 1.10"
    echo "juliaup default 1.10"
    echo "juliaup remove release"
    exit 1
fi

# Check if Julia version is 1.10
if ! julia --version | grep -q "1.10"; then
    echo "Julia 1.10 is not installed. Current version is $(julia --version)."
    echo "Please install Julia 1.10 by running:"
    echo "juliaup add 1.10"
    echo "juliaup default 1.10"
    echo "juliaup remove release"
    exit 1
fi

# Print the path of the Julia binary
julia_path=$(which julia)
echo "Julia 1.10 binary path: $julia_path"

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
