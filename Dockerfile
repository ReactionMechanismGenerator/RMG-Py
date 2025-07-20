# Parent Image
FROM condaforge/miniforge3:latest

# Install Bash shell
RUN ln -snf /bin/bash /bin/sh

# Install system dependencies
#
# List of deps and why they are needed:
#  - make, gcc, g++ for building RMG
#  - git for downloading RMG respoitories
#  - wget for downloading conda install script
#  - libxrender1 required by RDKit
#  - ca-certificates added for HTTPS downloads
RUN apt-get update && \
    apt-get install -y \
        make \
        gcc \
        wget \
        git \
        g++ \
        libxrender1 \
        ca-certificates && \
    apt-get autoremove -y && \
    apt-get clean -y


# Install Julia 1.10 using juliaup
RUN wget -qO- https://install.julialang.org | sh -s -- --yes --default-channel 1.10 && \
    /root/.juliaup/bin/juliaup add 1.10 && \
    /root/.juliaup/bin/juliaup default 1.10 && \
    /root/.juliaup/bin/juliaup list && \
    rm -rf /root/.juliaup/downloads /root/.juliaup/tmp
ENV PATH="/root/.juliaup/bin:$PATH"

# Set Bash as the default shell for following commands
SHELL ["/bin/bash", "-c"]

# Add build arguments for RMG-Py, RMG-database, and RMS branches.
# The defaults are set here, but they can be overridden at build time
# using the --build-arg option, or in the continous integration CI.yml file.
ARG RMG_Py_Branch=main
ARG RMG_Database_Branch=main
ARG RMS_Branch=for_rmg

# cd
WORKDIR /rmg

# Clone the RMG base and database repositories
RUN git clone --single-branch --branch ${RMG_Py_Branch} --depth 1 https://github.com/ReactionMechanismGenerator/RMG-Py.git && \
    git clone --single-branch --branch ${RMG_Database_Branch} --depth 1 https://github.com/ReactionMechanismGenerator/RMG-database.git

WORKDIR /rmg/RMG-Py

# build the conda environment
RUN conda env create --file environment.yml
# Remove conda package cache to reduce image size
RUN rm -rf /miniconda/pkgs

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

# Build RMG
RUN make

# Install and link Julia dependencies for RMS
# setting this env variable fixes an issue with Julia precompilation on Windows
ENV JULIA_CPU_TARGET="x86-64,haswell,skylake,broadwell,znver1,znver2,znver3,cascadelake,icelake-client,cooperlake,generic"
ENV RMS_BRANCH=${RMS_Branch}
# Usually this is set automatically, but we're not actually running
# in an active conda environment when building the Docker so we need to set it manually
ENV PYTHON_JULIAPKG_PROJECT="/miniconda/envs/rmg_env/julia_env"
RUN source install_rms.sh

# RMG-Py should now be installed and ready - trigger precompilation and test run
RUN python rmg.py examples/rmg/rms_constant_V/input.py
# delete the results, restore input.py from git
RUN rm -rf examples/rmg/rms_constant_V/* && \
    git checkout -- examples/rmg/rms_constant_V/

# when running this image, open an interactive bash terminal inside the rmg_env conda environment
RUN sed -i 's/conda activate base/conda activate rmg_env/' ~/.bashrc
ENTRYPOINT ["/bin/bash", "--login"]

