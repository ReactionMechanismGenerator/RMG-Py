# Parent Image
FROM ubuntu:latest

# Install Bash shell
RUN ln -snf /bin/bash /bin/sh

# Install system dependencies
#
# List of deps and why they are needed:
#  - make, gcc, g++ for building RMG
#  - git for downloading RMG respoitories
#  - wget for downloading conda install script
#  - libxrender1 required by RDKit
RUN apt-get update
RUN apt-get install -y make gcc wget git g++ libxrender1

# Install conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda 
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="$PATH:/miniconda/bin"

# Set solver backend to mamba for speed
RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba

# Set Bash as the default shell for following commands
SHELL ["/bin/bash", "-c"]

# cd
WORKDIR /rmg

# Clone the RMG base and database repositories
RUN git clone -b main https://github.com/ReactionMechanismGenerator/RMG-Py.git
RUN git clone -b main https://github.com/ReactionMechanismGenerator/RMG-database.git

# build the conda environment
WORKDIR /rmg/RMG-Py
RUN echo -e "\n  - libstdcxx-ng<13\n" >> environment.yml  # patch for linux
RUN conda env create --file environment.yml

# This runs all subsequent commands inside the rmg_env conda environment
#
# Analogous to just activating the environment, which we can't actually do here
# since that requires running conda init and restarting the shell (not possible
# in a Dockerfile build script)
SHELL ["conda", "run", "--no-capture-output", "-n", "rmg_env", "/bin/bash", "-c"]

# Set environment variables as directed in the RMG installation instructions
ENV RUNNER_CWD=/rmg
ENV PYTHONPATH="$RUNNER_CWD/RMG-Py:$PYTHONPATH"
ENV PATH="$RUNNER_CWD/RMG-Py:$PATH"

# Build RMG
RUN make

# Install and link Julia dependencies for RMS
RUN python -c "import julia; julia.install(); import diffeqpy; diffeqpy.install()" || true
RUN julia -e 'using Pkg; Pkg.add(PackageSpec(name="ReactionMechanismSimulator",rev="main")); using ReactionMechanismSimulator' || true 

# RMG-Py should now be installed and ready - trigger precompilation and test run
RUN python-jl rmg.py examples/rmg/minimal/input.py
# delete the results, preserve input.py
RUN mv examples/rmg/minimal/input.py .
RUN rm -rf examples/rmg/minimal/*
RUN mv input.py examples/rmg/minimal/

# when running this image, open an interactive bash terminal inside the conda environment
RUN echo "source activate rmg_env" >> ~/.bashrc

# add a warning so that users will not try and use latest
# comment/remove this line when building release images
RUN echo "NOTE: This docker image is intended for developers and is NOT guaranteed to work, please pull the most recent official release of rmg (3.1.1) instead." >> ~/.bashrc

ENTRYPOINT ["/bin/bash", "--login"]
