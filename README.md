# <img align="top" src="https://raw.githubusercontent.com/ReactionMechanismGenerator/RMG-Py/main/documentation/source/_static/rmg-logo-small.png"> Reaction Mechanism Generator (RMG)

[![Build status](https://img.shields.io/travis/ReactionMechanismGenerator/RMG-Py/main.svg)](https://travis-ci.org/ReactionMechanismGenerator/RMG-Py)
[![Codecov report](https://img.shields.io/codecov/c/github/ReactionMechanismGenerator/RMG-Py/main.svg)](https://codecov.io/gh/ReactionMechanismGenerator/RMG-Py)
[![Codacy report](https://img.shields.io/codacy/grade/5c12cecf3d01400a92ea20e14ca0b880/main.svg)](https://www.codacy.com/app/ReactionMechanismGenerator/RMG-Py/dashboard)
[![GitHub release](https://img.shields.io/github/release/ReactionMechanismGenerator/RMG-Py.svg)](https://github.com/ReactionMechanismGenerator/RMG-Py/releases)
[![Anconda](https://img.shields.io/conda/v/rmg/rmg.svg)](https://anaconda.org/rmg/rmg)
[![Gitter](https://img.shields.io/gitter/room/ReactionMechanismGenerator/RMG-Py.svg)](https://gitter.im/ReactionMechanismGenerator/RMG-Py)
[![RMG Website](https://img.shields.io/website-up-down-green-red/http/rmg.mit.edu.svg?label=rmg%20website)](https://rmg.mit.edu/)

## Description
This repository contains the Python version of **Reaction Mechanism Generator (RMG)**,
a tool for automatically generating chemical reaction
mechanisms for modeling reaction systems including pyrolysis, combustion,
atmospheric science, and more.

It also includes **Arkane**, the package for calculating thermodynamics, high-pressure-limit
rate coefficients, and pressure dependent rate coefficients from quantum chemical calculations.
Arkane is compatible with a variety of ab initio quantum chemistry software programs:
Gaussian, Q-Chem, Molpro, Orca, Psi4, and TeraChem.

## Source Code Repository
- [RMG Github Repository](https://github.com/ReactionMechanismGenerator/RMG-Py): contains the latest source code for RMG
- [RMG-database Github Repository](https://github.com/ReactionMechanismGenerator/RMG-database): contains source code for the latest version of the database

## How to Install
You can either download the source from GitHub and compile yourself, or download the binaries from Anaconda.
Please see the [Download and Install](http://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/installation/index.html) page for detailed instructions.

## Documentation
- [RMG Documentation](http://ReactionMechanismGenerator.github.io/RMG-Py/users/rmg/index.html) ([PDF version](https://github.com/ReactionMechanismGenerator/RMG-Py/raw/main/documentation/RMG-Py_and_Arkane_Documentation.pdf))
- [Arkane Documentation](http://ReactionMechanismGenerator.github.io/RMG-Py/users/arkane/index.html) ([PDF version](https://github.com/ReactionMechanismGenerator/RMG-Py/raw/main/documentation/RMG-Py_and_Arkane_Documentation.pdf))
- [RMG API Reference](http://reactionmechanismgenerator.github.io/RMG-Py/reference/index.html) ([PDF version](https://github.com/ReactionMechanismGenerator/RMG-Py/raw/main/documentation/RMG-Py_API_Reference.pdf))

## How to Contribute
Please see the [Contributor Guidelines](https://github.com/ReactionMechanismGenerator/RMG-Py/wiki/RMG-Contributor-Guidelines)
for details on how to contribute to RMG-Py, Arkane, or RMG-database.

## How to Give Feedback

Please post any issues you may have to the [issues page](https://github.com/ReactionMechanismGenerator/RMG-Py/issues/)
or drop in to the [chat room](https://gitter.im/ReactionMechanismGenerator/RMG-Py) or email [rmg_dev@mit.edu](mailto:rmg_dev@mit.edu) if you have questions.  

## Useful Links
- [Interactive Website](https://rmg.mit.edu): Visit this site to visualize RMG-generated models, view the databases, and 
perform thermodynamics and kinetics searches
- [Wiki](https://github.com/ReactionMechanismGenerator/RMG-Py/wiki): a wiki for developer notes
- [Issues Page](https://github.com/ReactionMechanismGenerator/RMG-Py/issues/): view current issues and feature requests

## Credits
- [Professor William H. Green's research group](http://cheme.scripts.mit.edu/green-group/) at the 
[Massachusetts Institute of Technology](http://web.mit.edu/) 
- [Professor Richard H. West's research group](http://www.northeastern.edu/comocheng/) at 
[Northeastern University](http://www.northeastern.edu/). 

## Resources and References
The resources and relevant publications are listed [here](https://rmg.mit.edu/resources) on the RMG-website. 
Please at least cite our latest publication on Reaction Mechanism Generator v3.0 and other
relevant publications when publishing the results using our software.

## How to cite
Please include the following citations if RMG, RMG-database, and/or Arkane were used for an academic study:
- C.W. Gao, J.W. Allen, W.H. Green, R.H. West,
  [Reaction Mechanism Generator: Automatic construction of chemical kinetic mechanisms](https://doi.org/10.1016/j.cpc.2016.02.013),
  Computer Physics Communications 2016, 203, 212-225.
- M. Liu, A. Grinberg Dana, M.S. Johnson, M.J. Goldman, A. Jocher, A.M. Payne, C.A. Grambow, K. Han, N.W. Yee,
  E.J. Mazeau, K. Blondal, R.H. West, C.F. Goldsmith, W.H. Green,
  [Reaction Mechanism Generator v3.0: Advances in Automatic Mechanism Generation](https://doi.org/10.1021/acs.jcim.0c01480),
  Journal of Chemical Information and Modeling 2021, 61(6), 2686-2696.
- M. S. Johnson, X. Dong, A. Grinberg Dana, Y. Chung, D. Farina, R. J. Gillis, M. Liu, N. W. Yee, K. Blondal, 
  E. Mazeau, C. A. Grambow, A. M. Payne, K. A. Spiekermann, H.-W. Pang, C. F. Goldsmith, R. H. West, W. H. Green,
  [RMG Database for Chemical Property Prediction](https://pubs.acs.org/doi/10.1021/acs.jcim.2c00965),
  Journal of Chemical Information and Modeling 2022, 62(20), 4906–4915.

## License Information
RMG is a free, open-source software package (distributed under the [MIT/X11 license](https://github.com/ReactionMechanismGenerator/RMG-Py/blob/main/LICENSE.txt)).
