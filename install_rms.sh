# install_rms.sh
#
# Convenience script to install ReactionMechanismSimulator into an rmg_env conda environment
#
# Usage:
#   source install_rms.sh
# Options:
#   RMS_INSTALLER environment variable can be set to "continuous", "standard", or "developer"
#     "continuous": non-interactive install for CI environments
#     "standard": interactive install with user confirmation (default)
#     "developer": install from local RMS source code (requires RMS_PATH to be set)
#   RMS_PATH environment variable (required if RMS_INSTALLER="developer")
#     Path to local ReactionMechanismSimulator.jl source code for developer mode
#   RMS_BRANCH environment variable (optional), for "standard" or "continuous" modes
#     Git branch of ReactionMechanismSimulator.jl to install (default: for_rmg)
# Example:
#   RMS_INSTALLER=developer RMS_PATH=/path/to/ReactionMechanismSimulator.jl source install_rms.sh


# Defaults to "standard" if not already set via RMS_INSTALLER env variable
RMS_INSTALLER=${RMS_INSTALLER:-standard}

# RMS branch for standard or continuous installs. Set to "for_rmg" by default.
RMS_BRANCH=${RMS_BRANCH:-for_rmg}

# Get local RMS path if in developer mode
if [ "$RMS_INSTALLER" = "developer" ]; then
    echo "Using developer mode for RMS installation"
    # Check if RMS_PATH is set
    if [ -z "$RMS_PATH" ]; then
        printf "Please enter full path to your local RMS source code: "
        read -r RMS_PATH
    fi
    # Validate Project.toml exists
    if [ ! -f "$RMS_PATH/Project.toml" ]; then
        echo "ERROR: '$RMS_PATH' does not contain a Project.toml file."
        echo "Please set RMS_PATH to a valid ReactionMechanismSimulator.jl directory."
        return 1
    fi
    echo "Using local RMS path: $RMS_PATH"
fi

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

# Ask the user to confirm RMS is being installed in the correct conda environemnt.
# Skip confirmation if this is run under continuous mode.
if [ "$RMS_INSTALLER" = "continuous" ]; then
    echo "Current conda environment: $current_env"
else
    echo "    Please confirm that you want to install RMS into the current conda environment: '$current_env'"
    echo "    If this is not correct, abort and activate the correct environment before rerunning."
    printf "Proceed with installation in '%s'? (y/N): " "$current_env"
    read confirm
    case "$confirm" in
    [yY][eE][sS]|[yY])
        echo "✅ Proceeding with installation in '$current_env'"
        ;;
    *)
        echo "❌ Aborted. Please activate the correct conda environment and try again."
        return 1
        ;;
    esac
fi

# Set environment variables for the current environment, for future uses
# https://juliapy.github.io/PythonCall.jl/stable/pythoncall/#If-you-already-have-Python-and-required-Python-packages-installed
# Find the actual conda executable (not the shell function)
CONDA_EXE="${CONDA_EXE:-$(conda info --base)/bin/conda}"
conda env config vars set JULIA_CONDAPKG_BACKEND=Null
conda env config vars set JULIA_CONDAPKG_EXE="$CONDA_EXE"
conda env config vars set JULIA_PYTHONCALL_EXE=$CONDA_PREFIX/bin/python
conda env config vars set PYTHON_JULIAPKG_EXE=$(which julia)
conda env config vars set PYTHON_JULIAPKG_PROJECT=$CONDA_PREFIX/julia_env
# Also export for current shell/session (needed for Docker/non-interactive use)
export JULIA_CONDAPKG_BACKEND=Null
export JULIA_CONDAPKG_EXE="$CONDA_EXE"
export JULIA_PYTHONCALL_EXE="$CONDA_PREFIX/bin/python"
export PYTHON_JULIAPKG_EXE="$(which julia)"
export PYTHON_JULIAPKG_PROJECT="$CONDA_PREFIX/julia_env"

echo "Environment variables referencing JULIA:"
env | grep JULIA

# Initialize the Julia environment from Python using juliacall
python << EOF || return 1
import sys
try:
    from juliacall import Main
    Main.seval('println("Active Julia environment: ", Base.active_project())')
    Main.seval('println("Julia load path: ", Base.load_path())')
    Main.seval('using Pkg')
    Main.seval('Pkg.status()')
except Exception as e:
    print("❌ Error while initializing Julia environment:")
    print(e)
    sys.exit(1)
EOF

# Install RMS
if [ "$RMS_INSTALLER" = "standard" ] || [ "$RMS_INSTALLER" = "continuous" ]; then
    echo "Installing RMS from branch: $RMS_BRANCH"
    julia << 'EOF' || echo "RMS standard install error - continuing anyway ¯\\_(ツ)_/¯"
    using Pkg
    Pkg.activate(ENV["PYTHON_JULIAPKG_PROJECT"])
    Pkg.add(Pkg.PackageSpec(name="ReactionMechanismSimulator", url="https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl.git", rev=ENV["RMS_BRANCH"]))
    println(read(joinpath(dirname(pathof(ReactionMechanismSimulator)), "..", "Project.toml"), String))
    Pkg.instantiate()
    try
        @info "Loading RMS"
        using ReactionMechanismSimulator
        @info "RMS loaded successfully!"
    catch err
        @error "Failed to load RMS" exception=err
        Base.show_backtrace(stderr, catch_backtrace())
        exit(1)
    end
EOF
elif [ "$RMS_INSTALLER" = "developer" ]; then
    echo "Installing RMS in developer mode from path: $RMS_PATH"
    julia << 'EOF' || echo "RMS developer install error - continuing anyway ¯\\_(ツ)_/¯"
    using Pkg
    println(ENV["PYTHON_JULIAPKG_PROJECT"])
    Pkg.activate(ENV["PYTHON_JULIAPKG_PROJECT"])
    Pkg.develop(path=ENV["RMS_PATH"])
    Pkg.instantiate()
    try
        @info "Loading RMS"
        using ReactionMechanismSimulator
        @info "RMS loaded successfully!"
    catch err
        @error "Failed to load RMS" exception=err
        Base.show_backtrace(stderr, catch_backtrace())
        exit(1)
    end
EOF
else
    echo "Unknown RMS_INSTALLER mode: $RMS_INSTALLER. Must be either 'continuous', 'standard' or 'developer'."
    return 1
fi

julia_status=$?
if [ $julia_status -ne 0 ]; then
    echo "RMS installation failed!"
    return $julia_status
fi

# this makes the above RMS installation available to Python
conda install -y conda-forge::pyjuliacall

echo "Checking if ReactionMechanismSimulator is installed in the current conda environment for Python usage..."

python << EOF
import sys
try:
    from juliacall import Main
    RMS_Pkg = Main.seval('Base.identify_package("ReactionMechanismSimulator")')
    print("Package identify result: ", RMS_Pkg)
    if RMS_Pkg is Main.nothing:
        print("❌ ReactionMechanismSimulator is NOT installed correctly.")
        sys.exit(1)
    else:
        print("✅ ReactionMechanismSimulator is succesfully installed!")
        sys.exit(0)
except Exception as e:
    print("❌ Error while checking ReactionMechanismSimulator installation:")
    print(e)
    sys.exit(1)
EOF
