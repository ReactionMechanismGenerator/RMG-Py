#!/bin/bash

# Set up anaconda
wget http://repo.continuum.io/miniconda/Miniconda-2.2.2-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b
export PATH=/home/travis/anaconda/bin:$PATH

# Update conda itself
conda update --yes conda