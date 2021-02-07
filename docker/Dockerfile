# Set base image to miniconda to save space
FROM continuumio/miniconda3:4.9.2

# Update conda
RUN conda update -n base -c defaults conda

# Install make and gcc so we can compile RMG
RUN apt-get update && apt-get install -y make gcc g++ git

# Create RMG user within image
RUN useradd --create-home rmg

# Switch to RMG user
USER rmg

# Add conda setup to .bashrc for interactive use (does not affect Docker build)
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/rmg/.bashrc && \
    echo "conda activate rmg_env" >> /home/rmg/.bashrc

# Set up PATH (persistent in final image)
ENV PATH /home/rmg/.conda/envs/rmg_env/bin:/opt/conda/bin:/home/rmg/RMG-Py:$PATH

# Set working directory to /rmg
WORKDIR /home/rmg

# Get environment file from GitHub and create environment
# Do this here to take better advantage of build caching since the environment does not change as often as RMG
ADD --chown=rmg:rmg https://raw.githubusercontent.com/ReactionMechanismGenerator/RMG-Py/stable/environment.yml .
RUN conda env create -f environment.yml

# Clone repos
RUN git clone https://github.com/ReactionMechanismGenerator/RMG-Py.git && \
    git clone https://github.com/ReactionMechanismGenerator/RMG-database.git

# Compile RMG
RUN cd RMG-Py && make

# Set default command as rmg.py --help
CMD ["rmg.py", "--help"]

