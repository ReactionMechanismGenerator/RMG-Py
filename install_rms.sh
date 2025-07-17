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
current_env=$(conda info --envs | grep -v '^#' | awk '/\*/{print $1}')
echo "Current conda environment: $current_env"

# Set environment variables for the current environment, for future uses
# https://juliapy.github.io/PythonCall.jl/stable/pythoncall/#If-you-already-have-Python-and-required-Python-packages-installed
conda env config vars set JULIA_CONDAPKG_BACKEND=Null
conda env config vars set JULIA_CONDAPKG_EXE=$(which conda)
conda env config vars set JULIA_PYTHONCALL_EXE=$CONDA_PREFIX/bin/python
conda env config vars set PYTHON_JULIAPKG_EXE=$(which julia)
conda env config vars set PYTHON_JULIAPKG_PROJECT=$CONDA_PREFIX/julia_env
# Also export for current shell/session (needed for Docker/non-interactive use)
export JULIA_CONDAPKG_BACKEND=Null
export JULIA_CONDAPKG_EXE="$(which conda)"
export JULIA_PYTHONCALL_EXE="$CONDA_PREFIX/bin/python"
export PYTHON_JULIAPKG_EXE="$(which julia)"
export PYTHON_JULIAPKG_PROJECT="$CONDA_PREFIX/julia_env"

conda install -y conda-forge::pyjuliacall

echo "Environment variables referencing JULIA:"
env | grep JULIA

# Ask user whether to do a standard or developer install
if [ -z "$RMS_MODE" ]; then
    echo "Choose installation mode:"
    echo "  1) Standard install (download from GitHub)"
    echo "  2) Developer install (install from local path)"
    read -p "Enter 1 or 2 [default: 1]: " installation_choice

    if [ "$installation_choice" = "2" ]; then
        RMS_MODE="dev"
    else
        RMS_MODE="standard"
    fi
fi

echo "Selected RMS installation mode: $RMS_MODE"

# Default RMS branch for standard install
RMS_BRANCH=${RMS_BRANCH:-for_rmg}

# Ask for local RMS path
if [ "$RMS_MODE" = "dev" ]; then
    read -e -p "Please enter full path to your local RMS source code: " RMS_PATH
    if [ ! -d "$RMS_PATH" ]; then
        echo "ERROR: '$RMS_PATH' is not a valid directory."
        exit 1
    fi
    echo "Using local RMS path: $RMS_PATH"
fi

if [ "$RMS_MODE" = "standard" ]; then
    echo "Installing RMS from branch: $RMS_BRANCH"
    julia << EOF || echo "RMS standard install error - continuing anyway ¯\\_(ツ)_/¯"
    using Pkg
    Pkg.activate(ENV["PYTHON_JULIAPKG_PROJECT"])
    Pkg.add(Pkg.PackageSpec(name="ReactionMechanismSimulator", url="https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl.git", rev="$RMS_BRANCH"))
    println(read(joinpath(dirname(pathof(ReactionMechanismSimulator)), \"..\", \"Project.toml\"), String))
    Pkg.instantiate()
    using ReactionMechanismSimulator
EOF
elif [ "$RMS_MODE" = "dev" ]; then
    echo "Installing RMS in developer mode from path: $RMS_PATH"
    julia << EOF || echo "RMS developer install error - continuing anyway ¯\\_(ツ)_/¯"
    using Pkg
    Pkg.activate(ENV["PYTHON_JULIAPKG_PROJECT"])
    Pkg.develop(path="$RMS_PATH")
    Pkg.instantiate()
    using ReactionMechanismSimulator
EOF
else
    echo "Unknown RMS_MODE: $RMS_MODE. Must be either 'standard' or 'dev'."
    exit 1
fi

echo "Checking if ReactionMechanismSimulator is installed in the current conda environment for Python usage..."
python -c "from juliacall import Main; import sys; sys.exit(0 if Main.seval('Base.identify_package(\"ReactionMechanismSimulator\") !== nothing') and print('ReactionMechanismSimulator is installed in $current_env') is None else 1)"
