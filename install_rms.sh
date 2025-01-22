# manual installation steps
# start by installing juliaup on your system globally
# curl -fsSL https://install.julialang.org | sh
# # restart shell
# juliaup add 1.9
# juliaup default 1.9
# juliaup remove release

# actual steps
conda install -y conda-forge::pyjuliacall  # conda-forge::pyside2
# export JULIA_CONDAPKG_BACKEND=Current
# export JULIA_CONDAPKG_EXE=$CONDA_EXE
# export JULIA_PYTHONCALL_EXE=$CONDA_PREFIX/bin/python

# https://juliapy.github.io/PythonCall.jl/stable/pythoncall/#If-you-already-have-Python-and-required-Python-packages-installed
export JULIA_CONDAPKG_BACKEND=Null
export JULIA_PYTHONCALL_EXE=$CONDA_PREFIX/bin/python
# export MPLBACKEND=qt5agg  # supported backend for PythonPlot.jl, needed by RMS

julia -e 'using Pkg; Pkg.add(Pkg.PackageSpec(name="ReactionMechanismSimulator", url="https://github.com/hwpang/ReactionMechanismSimulator.jl.git", rev="fix_installation")); using ReactionMechanismSimulator; Pkg.instantiate()'

# ensure that juliacall in Python uses the correct julia executable: https://github.com/JuliaPy/PyJuliaPkg?tab=readme-ov-file#which-julia-gets-used
export PYTHON_JULIAPKG_EXE=$(which julia)
