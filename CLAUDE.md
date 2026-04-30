# CLAUDE.md

Guidance for AI agents working in this repo. Also read [.github/copilot-instructions.md](.github/copilot-instructions.md) for project overview, package layout, and code patterns. This file focuses on the non-obvious mechanics of building, testing, and contributing.

## What this repo is

**RMG-Py** = Reaction Mechanism Generator (chemical kinetics mechanism generator). Two installable packages: `rmgpy/` (mechanism generation) and `arkane/` (statistical mechanics from QM calculations). Heavy use of Cython for performance. Depends on a sibling repo **RMG-database** for thermo/kinetics data.

Entry points:
- `rmg.py <input.py>` — runs RMG on a Python-style input file (camelCase syntax — see `examples/rmg/minimal/input.py`)
- `Arkane.py <input.py>` — runs Arkane
- Both are also installed as console scripts (`rmg.py`, `Arkane.py`) via `setup.py`. They thin-wrap `rmgpy.__main__:main` / `arkane.__main__:main`.
- Scripts under `scripts/` are also installed (e.g. `simulate.py`, `diffModels.py`, `mergeModels.py`, `rmg2to3.py`).

## Setup

The conda env is the only supported install path. Python is pinned `>=3.9,<3.12`.

```bash
conda env create --file environment.yml   # creates env named `rmg_env`
conda activate rmg_env
git clone https://github.com/ReactionMechanismGenerator/RMG-database ../RMG-database
make install
```

`make install` runs `python utilities.py check-pydas` (which writes [rmgpy/solver/settings.pxi](rmgpy/solver/settings.pxi) — see Cython section), then `pip install --no-build-isolation -vv -e .`, then touches a `.installed` sentinel. Subsequent `make` invocations skip reinstall unless the sentinel is missing.

**Always keep [environment.yml](environment.yml) and [.conda/meta.yaml](.conda/meta.yaml) in sync** — both define runtime deps and CI builds from `meta.yaml` for the conda package.

Optional pieces:
- `./install_rms.sh` — installs ReactionMechanismSimulator (Julia-based reactor backend). Required for `rms*` reactor types in input files. Honors `RMS_INSTALLER={continuous,standard,developer}` and `RMS_BRANCH` (default `for_rmg`).
- `make q2dtor` — clones Q2DTor into `external/` for 2D rotor calculations in Arkane.

## Build / Cython

Cython modules are listed explicitly in `setup.py` `ext_modules`. **Some `.py` files are cythonized** (not just `.pyx`): e.g. `rmgpy/molecule/molecule.py`, `group.py`, `atomtype.py`, `rmgpy/species.py`, `rmgpy/reaction.py`, `rmgpy/quantity.py`, `rmgpy/constants.py`. If you edit one of these, **rebuild** — the `.so` is what gets imported, not the `.py`.

Workflow:
- `make build` — incremental in-place `setup.py build_ext --inplace`. Fast. Use this after editing `.pyx`/`.pxd`/cythonized `.py`.
- `make` (default `all`) — checks deps, ensures `.installed` sentinel, then `make build`. Safe go-to.
- `make clean` — removes `.so`, `.pyc`, generated `.c`, `build/`, and `.installed`. Also `pip uninstall`s the package.
- `make decython` — deletes most `.so` files (keeps `_statmech.so`, `quantity.so`, and `rmgpy/solver/*.so`) so pure Python is loaded for debugging. **Pure Python mode is not reliably tested**; expect breakage.

Cython conventions in this repo:
- Compile language level is Python 3.
- Pair every public `cdef class` / `cpdef` method with a `.pxd` declaration.
- New extension files **must be added to `ext_modules` in `setup.py`** or they will silently not be built.
- The DASPK/DASSL solver is selected at compile time via `rmgpy/solver/settings.pxi` (auto-written by `utilities.py check-pydas` from whatever PyDAS variant is installed). Do not commit changes to `settings.pxi`.
- macOS-specific: `setup.py` deduplicates `-Wl,-rpath` flags from sysconfig before invoking Cython, to work around an LC_RPATH issue with conda-forge's Python on darwin. Don't remove that block.

## Tests

Configured in [pytest.ini](pytest.ini): `testpaths = test`, `python_files = *Test.py`, `python_classes = *Test Test*`. Tests live under `test/` mirroring `rmgpy/` and `arkane/`.

Default pytest flags include `-s -vv --keep-duplicates` and coverage (`--cov=arkane --cov=rmgpy --cov-report html`). `test/regression/` is excluded.

Markers (from `pytest.ini`):
- `@pytest.mark.functional` — slower functional tests
- `@pytest.mark.database` — tests that require RMG-database to be cloned and loaded
- Unmarked = unit tests

Make targets:
```bash
make test            # unit tests only (excludes functional, database)
make test-functional
make test-database
make test-all        # everything
```

Run a subset directly:
```bash
pytest test/rmgpy/molecule/atomtypeTest.py
pytest -k "test_pattern"
pytest -m "functional"
```

`pytest-xdist` (`-n auto`) is supported but **incompatible with RMS/Julia** — only use when RMS is not installed.

`test/conftest.py` forces `multiprocessing.set_start_method('fork')` and silences OpenBabel error logging. Be aware of the `fork` start method when adding tests that touch multiprocessing.

### Regression tests

Separate from pytest. Each `test/regression/<name>/` has an `input.py`. CI runs `python rmg.py test/regression/<name>/input.py` and diffs core/edge models against artifacts produced on `main`. Locally you can reproduce a single one:
```bash
python rmg.py test/regression/superminimal/input.py
python scripts/checkModels.py ...   # (see .github/workflows/CI.yml for arg shape)
```
Adding a new regression test means editing the **two lists** in [.github/workflows/CI.yml](.github/workflows/CI.yml) (Execution + Comparison steps); the first PR will fail CI until baseline artifacts exist on `main`.

The `Makefile` also has `eg0`-`eg10` targets that copy example inputs into `testing/<name>/` and run `rmg.py` — useful for ad-hoc end-to-end smoke testing (`eg0` is fastest).

## Linting / formatting / typing

There is **no configured linter, formatter, or type checker** in this repo (no `pyproject.toml`, `ruff.toml`, `.flake8`, `mypy.ini`, or `pre-commit` config). The only style guidance is "follow PEP 8 for new code, but don't churn existing code just for style." Don't run `black`/`ruff format`/`isort` over the tree as part of unrelated changes — diffs balloon and reviews stall.

## Database integration

RMG looks up `database.directory` in this order:
1. `database.load(path=...)` arg in code
2. `rmgrc` in cwd
3. `~/.rmg/rmgrc`
4. `rmgpy/rmgrc` (alongside the package)
5. Default: `../RMG-database/input` relative to RMG-Py source

Template: [rmgpy/rmgrc_template](rmgpy/rmgrc_template). Copy it (don't edit in place — it's overwritten on install). In CI, the database is checked out at the branch named in `RMG_DATABASE_BRANCH` (env var in [.github/workflows/CI.yml](.github/workflows/CI.yml)); change that line if your PR depends on an unmerged database branch.

## Conventions

- **Python API uses `snake_case`**. The mass rename happened in the Python 3 transition (see `scripts/rmg2to3.py` — automated converter for old code).
- **Input file DSL keeps `camelCase`** (`thermoLibraries`, `simpleReactor`, `terminationConversion`, ...) for backward compatibility. When adding a new input keyword, follow camelCase and update [documentation/source/users/rmg/input.rst](documentation/source/users/rmg/input.rst).
- All source files require the MIT license header (template lives in [LICENSE.txt](LICENSE.txt); `python utilities.py update-headers` re-applies it across `.py`/`.pyx`/`.pxd` in `rmgpy/`, `scripts/`, and the root).
- Use `logging` not `print`.
- Don't reach for `__init__.py`-as-namespace imports across cython modules; use `cimport rmgpy.constants as constants` etc.

## Documentation

Sphinx docs in `documentation/source/`. Build with `make documentation` (calls `make -C documentation html`). Output: `documentation/build/html/index.html`. Built with the standard `rmg_env` (no separate doc env).

When changing things, also update:
- **Input file syntax / new options** → [documentation/source/users/rmg/input.rst](documentation/source/users/rmg/input.rst) (this is treated as required by reviewers).
- **New public API** → ensure docstrings exist; add module to a toctree under `documentation/source/reference/` if the module is new. API docs are auto-generated via `sphinx.ext.autodoc`.
- **New user-facing feature** → mention in `documentation/source/users/rmg/features.rst` or a sibling `.rst`.
- **Behavior change** → relevant section of the user guide (`users/rmg/` or `users/arkane/`).

The `gh-pages` branch hosts the live site; CI publishes on push to `main`.

## CI

- [.github/workflows/CI.yml](.github/workflows/CI.yml): build + tests on Linux (ubuntu-latest, all Python versions, with and without RMS) and macOS (latest Python only). Linux runs `make test-all`; other matrix entries run only unit tests. Regression job runs separately on `ubuntu-latest`.
- [.github/workflows/conda_build.yml](.github/workflows/conda_build.yml): builds the conda package from [.conda/meta.yaml](.conda/meta.yaml).
- [.github/workflows/docs.yml](.github/workflows/docs.yml): Sphinx build + publish to `gh-pages`.

## Quick gotchas

- **Edits to `.pyx`/`.pxd`/cythonized `.py` won't take effect until you rebuild** (`make build`). Mysterious unchanged behavior is almost always a stale `.so`.
- **`.so` files persist across branch switches.** When chasing a weird bug after a checkout, `make clean && make` before debugging.
- **Don't use `--no-verify` or skip Cython rebuilds** to make a commit go through; the underlying issue will resurface in CI.
- **Functional/database tests need RMG-database checked out** at a compatible branch in `../RMG-database`.
- **RMS reactor types in input files require Julia** — without `install_rms.sh` they'll fail at runtime, not import.
