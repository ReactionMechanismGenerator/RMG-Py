.. _installation:

************
Installation
************

.. NOTE::
    RMG has been tested on the Python 2.5, 2.6, and 2.7 releases; dependency issues render it incompatible with Python 3.x releases


Briefly, RMG depends on the following packages:

* **NumPy:** fast matrix operations
* **SciPy:** fast mathematical toolkit
* **matplotlib:** generating plots
* **rdkit:** open-source cheminformatics toolkit
* **guppy:** memory profiling tools
* **Cython:** compiling Python modules to C
* **quantities:** unit conversion
* **nose:** advanced unit test controls
* **Sphinx:** documentation generation
* **pydot:** interface to Dot graph language
* **cairo:** molecular diagram generation
* **psutil:** system utilization diagnostic tool
* **xlwt:** generating Excel output files
* **Graphviz:** generating flux diagrams
* **PyDAS:** differential algebraic system solver
* **PyDQED:** constrained nonlinear optimization
* **RMG-database:** thermodynamic and kinetic libraries

Refer to these platform-specific instructions for details on the best ways to install these packages before attempting to build RMG-Py:

.. toctree::
    :maxdepth: 1

    windows
    linux
    macos

See also the instructions for installing these optional components:

.. toctree::
    :maxdepth: 1

    QMthermo
