# manual installation steps
# start by installing juliaup on your system globally
# curl -fsSL https://install.julialang.org | sh
# # restart shell
# juliaup add 1.9
# juliaup default 1.9

# actual steps
conda install -y conda-forge::pyjuliacall
export JULIA_CONDAPKG_BACKEND=Current
export JULIA_CONDAPKG_EXE=$CONDA_EXE
export JULIA_PYTHONCALL_EXE=$CONDA_PREFIX/bin/python
julia -e 'using Pkg; Pkg.add(Pkg.PackageSpec(name="ReactionMechanismSimulator", url="https://github.com/hwpang/ReactionMechanismSimulator.jl.git", rev="fix_installation")); using ReactionMechanismSimulator'
