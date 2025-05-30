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
RUN apt-get update && \
    apt-get install -y \
        make \
        gcc \
        wget \
        git \
        g++ \
        libxrender1 && \
    apt-get autoremove -y && \
    apt-get clean -y

# Install conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="$PATH:/miniconda/bin"

# Set Bash as the default shell for following commands
SHELL ["/bin/bash", "-c"]

# Add build arguments for RMG-Py, RMG-database, and RMS branches.
ARG RMG_Py_Branch=main
ARG RMG_Database_Branch=main
ARG RMS_Branch=main
ENV rmsbranch=${RMS_Branch}

# cd
WORKDIR /rmg

# Clone the RMG base and database repositories
RUN git clone --single-branch --branch ${RMG_Py_Branch} --depth 1 https://github.com/ReactionMechanismGenerator/RMG-Py.git && \
    git clone --single-branch --branch ${RMG_Database_Branch} --depth 1 https://github.com/ReactionMechanismGenerator/RMG-database.git

WORKDIR /rmg/RMG-Py
# build the conda environment
RUN conda env create --file environment.yml

# This runs all subsequent commands inside the rmg_env conda environment
#
# Analogous to just activating the environment, which we can't actually do here
# since that requires running conda init and restarting the shell (not possible
# in a Dockerfile build script)
SHELL ["conda", "run", "--no-capture-output", "-n", "rmg_env", "/bin/bash", "-c"]

RUN conda clean --all --yes

# Set environment variables as directed in the RMG installation instructions
ENV RUNNER_CWD=/rmg
ENV PATH="$RUNNER_CWD/RMG-Py:$PATH"

# 1. Build RMG
# 2. Install and link Julia dependencies for RMS
# setting this env variable fixes an issue with Julia precompilation on Windows
ENV JULIA_CPU_TARGET="x86-64,haswell,skylake,broadwell,znver1,znver2,znver3,cascadelake,icelake-client,cooperlake,generic"
RUN source install_rms.sh

# RMG-Py should now be installed and ready - trigger precompilation and test run
RUN python rmg.py examples/rmg/minimal/input.py
# delete the results, preserve input.py
RUN mv examples/rmg/minimal/input.py . && \
    rm -rf examples/rmg/minimal/* && \
    mv input.py examples/rmg/minimal/

# when running this image, open an interactive bash terminal inside the conda environment
RUN echo "source activate rmg_env" >~/.bashrc
ENTRYPOINT ["/bin/bash", "--login"]
