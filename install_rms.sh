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
    return 1
fi

# Check if julia command is available
if ! command -v julia &> /dev/null; then
    echo "Julia is not installed. Please install version 1.10 by running:"
    echo "juliaup add 1.10"
    echo "juliaup default 1.10"
    echo "juliaup remove release"
    return 1
fi

# Check if Julia version is 1.10
if ! julia --version | grep -q " 1\.10"; then
    echo "Julia 1.10 is not installed. Current version is $(julia --version)."
    echo "Please install Julia 1.10 by running:"
    echo "juliaup add 1.10"
    echo "juliaup default 1.10"
    echo "juliaup remove release"
    return 1
fi

# Print the path of the Julia binary
julia_path=$(which julia)
echo "Julia 1.10 binary path: $julia_path"

# Get current conda environment name
current_env=$(conda info --envs | grep '\*' | awk '{print $1}')
echo "Current conda environment: $current_env"

# Set environment variables for the current environment, for future uses
# https://juliapy.github.io/PythonCall.jl/stable/pythoncall/#If-you-already-have-Python-and-required-Python-packages-installed
conda env config vars set JULIA_CONDAPKG_BACKEND=Null
conda env config vars set JULIA_PYTHONCALL_EXE=$CONDA_PREFIX/bin/python
conda env config vars set PYTHON_JULIAPKG_EXE=$(which julia)
conda env config vars set PYTHON_JULIAPKG_PROJECT=$CONDA_PREFIX/julia_env
# Also export for current shell/session (needed for Docker/non-interactive use)
export JULIA_CONDAPKG_BACKEND=Null
export JULIA_PYTHONCALL_EXE="$CONDA_PREFIX/bin/python"
export PYTHON_JULIAPKG_EXE="$(which julia)"
export PYTHON_JULIAPKG_PROJECT="$CONDA_PREFIX/julia_env"

conda install -y conda-forge::pyjuliacall

echo "Environment variables referencing JULIA:"
env | grep JULIA

# Use RMS_BRANCH environment variable if set, otherwise default to for_rmg
RMS_BRANCH=${RMS_BRANCH:-for_rmg}
echo "Installing ReactionMechanismSimulator from branch: $RMS_BRANCH"

julia -e "using Pkg; Pkg.add(Pkg.PackageSpec(name=\"ReactionMechanismSimulator\", url=\"https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl.git\", rev=\"$RMS_BRANCH\")); using ReactionMechanismSimulator; Pkg.instantiate()" || echo "RMS install error - continuing anyway ¯\_(ツ)_/¯"

echo "Checking if ReactionMechanismSimulator is installed in the current conda environment for Python usage..."
python -c "from juliacall import Main; import sys; sys.exit(0 if Main.seval('Base.identify_package(\"ReactionMechanismSimulator\") !== nothing') and print('ReactionMechanismSimulator is installed in $current_env') is None else 1)"
