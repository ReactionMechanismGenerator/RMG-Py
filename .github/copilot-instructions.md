# RMG-Py Copilot Instructions

## Project Overview
RMG-Py is the Reaction Mechanism Generator - an automatic chemical kinetics mechanism generator. It consists of two main components:
- **RMG** (`rmgpy/`): Core mechanism generation engine
- **Arkane** (`arkane/`): Statistical mechanics and transition state theory calculations

## Architecture

### Core Packages
- `rmgpy/molecule/` - Molecular graph representation (`Molecule`, `Atom`, `Bond`, `Group`)
- `rmgpy/thermo/` - Thermodynamic models (NASA, Wilhoit, ThermoData)
- `rmgpy/kinetics/` - Rate coefficient models (Arrhenius, Chebyshev, pressure-dependent)
- `rmgpy/solver/` - ODE solvers for reactor simulations
- `rmgpy/rmg/` - Main RMG algorithm (`main.py`, `model.py`, `react.py`)
- `rmgpy/data/` - Database interfaces for thermo, kinetics, transport

### Key Base Classes
- `RMGObject` (in `rmgpy/rmgobject.pyx`) - Base class providing `as_dict()`/`make_object()` for YAML serialization
- `Graph`/`Vertex`/`Edge` (in `rmgpy/molecule/graph.pyx`) - Graph isomorphism via VF2 algorithm
- `Species` and `Reaction` are central objects connecting molecules to thermodynamics and kinetics

### Cython Architecture
Performance-critical code uses Cython (`.pyx` files) with declaration files (`.pxd`):
- Always pair `.pyx` with `.pxd` for public cdef classes/methods
- Use `cpdef` for methods callable from both Python and Cython
- Use `cimport` for Cython-level imports (e.g., `cimport rmgpy.constants as constants`)
- Register new Cython modules in `setup.py` `ext_modules` list

## Development Commands
```bash
make install          # Build Cython extensions and install in editable mode
make test             # Run unit tests (excludes functional/database tests)
make test-functional  # Run functional tests
make test-database    # Run database tests
make test-all         # Run all tests
make clean            # Remove build artifacts
make decython         # Remove .so files for "pure Python" debugging. This target probably broken.
make documentation    # Build Sphinx docs
```

## Testing Conventions
- Tests live in `test/` mirroring `rmgpy/` and `arkane/` structure
- Test files: `*Test.py` (e.g., `speciesTest.py`, `reactionTest.py`)
- Test classes: `class TestClassName:` or `class ClassNameTest:`
- Use `pytest` with fixtures (`@pytest.fixture(autouse=True)` for setup)
- Markers: `@pytest.mark.functional`, `@pytest.mark.database`
- Run specific tests: `pytest -k "test_name_pattern"`

## Code Patterns

### Molecular Representations
```python
from rmgpy.molecule import Molecule
mol = Molecule().from_smiles("CC")  # From SMILES
mol = Molecule().from_adjacency_list("""...""")  # From adjacency list
mol.is_isomorphic(other_mol)  # Graph isomorphism check
```

### Species and Reactions
```python
from rmgpy.species import Species
species = Species(label='ethane', molecule=[Molecule().from_smiles("CC")])
species.generate_resonance_structures()
```

## Input Files
- RMG inputs: Python scripts defining `database()`, `species()`, `simpleReactor()`, etc.
- See `examples/rmg/minimal/input.py` for structure
- Arkane inputs: Python scripts with `species()`, `transitionState()`, `reaction()` blocks

## RMG-database Integration
The **RMG-database** is a separate repository containing all thermodynamic, kinetics, and transport data. It's typically cloned alongside RMG-Py in a sibling folder named `RMG-database`.

### Database Structure (in RMG-database repo)
- `thermo/` - Thermodynamic libraries and group additivity data
- `kinetics/families/` - Reaction family templates with rate rules (e.g., `H_Abstraction`, `R_Addition_MultipleBond`)
- `kinetics/libraries/` - Curated rate coefficient libraries
- `solvation/` - Solvent and solute parameters
- `transport/` - Transport properties

### How RMG-Py Loads the Database
The `RMGDatabase` class (`rmgpy/data/rmg.py`) is the central interface:
```python
from rmgpy.data.rmg import RMGDatabase
database = RMGDatabase()
database.load(
    path='/path/to/RMG-database',
    thermo_libraries=['primaryThermoLibrary'],
    kinetics_families='default',
    reaction_libraries=[],
)
```

### Key Database Classes
- `ThermoDatabase` (`rmgpy/data/thermo.py`) - Estimates thermo via group additivity or libraries
- `KineticsDatabase` (`rmgpy/data/kinetics/database.py`) - Manages reaction families and libraries
- `KineticsFamily` (`rmgpy/data/kinetics/family.py`) - Template-based reaction generation using `Group` pattern matching
- `Entry` (`rmgpy/data/base.py`) - Base class for database entries with metadata

### Data Flow for Species Thermodynamics
1. `Species.get_thermo_data()` → `ThermoDatabase.get_thermo_data(species)`
2. First checks thermo libraries for exact match (via graph isomorphism)
3. Falls back to group additivity estimation using functional group contributions
4. Returns `ThermoData`, `NASA`, or `Wilhoit` object

### Data Flow for Reaction Kinetics
1. `KineticsFamily.generate_reactions(reactants)` - Matches reactant molecules to family templates
2. Creates `TemplateReaction` objects with labeled atoms from template matching
3. `KineticsFamily.get_kinetics()` - Estimates rate using rate rules or training reactions
4. Returns `Arrhenius` or pressure-dependent kinetics model

## External Dependencies
- **RMG-database**: Set location via `RMG_DATABASE_BRANCH` env var in CI, or pass path to `database.load()`
- **Julia/RMS**: Optional (recommended) reactor simulation backend (install via `./install_rms.sh`)
- Environment managed via `environment.yml` (conda/mamba)

## Documentation
Documentation lives in `documentation/source/` and is built with Sphinx (`make documentation`).

### User Documentation (`documentation/source/users/`)
- `users/rmg/` - RMG user guide (how to run, configure, interpret output)
- `users/arkane/` - Arkane user guide
- **Critical file**: `users/rmg/input.rst` - Documents all input file options. **Must be updated when changing input file syntax or adding new features.**

### API Reference (`documentation/source/reference/`)
- Auto-generated from docstrings using `sphinx.ext.autodoc`
- Each module has a corresponding `.rst` file (e.g., `reference/species/index.rst` → `rmgpy/species.py`)
- **Maintenance**: Add new modules to the appropriate `index.rst` toctree. Docstrings in code are automatically extracted.
- Uses reStructuredText format with `.. automodule::` directives

### When to Update Documentation
- **New input file options**: Update `users/rmg/input.rst`
- **New public API**: Ensure docstrings exist; add module to `reference/` if new
- **Changed behavior**: Update relevant user guide section
- **New features**: Add to `users/rmg/features.rst` or create new `.rst` file

## Style Guidelines
- Follow PEP 8
- Docstrings describe purpose, not implementation
- Use `logging` module (not print statements)
- MIT/X11 license header required on all source files
