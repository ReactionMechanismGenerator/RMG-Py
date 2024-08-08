## How to contribute to RMG-Py

Thank you for contributing to RMG-Py! Please take a moment to review our guidelines:

### **Did you find a bug? Do you want to see a new feature?**

* Please open an Issue to the corresponding repository:
    * [RMG-Py](https://github.com/ReactionMechanismGenerator/RMG-Py/issues): For functionality of the RMG and Arkane software packages
    * [RMG-database](https://github.com/ReactionMechanismGenerator/RMG-database/issues): For issues related to the data available to RMG
    * [RMG-website](https://github.com/ReactionMechanismGenerator/RMG-website/issues): For issues related to the [RMG website](https://rmg.mit.edu)


### **Did you write code that fixes a bug or adds a new feature?**

* Open a new GitHub PR to merge into the main branch. Make sure the PR clearly describes the problem + solution. If applicable, include the relevant issue.

* Your PR must pass unit tests, regression tests, and code coverage, and receive approval from at least one reviewer before it can be merged in.

* If you wrote a new feature, please [add unit tests](https://github.com/ReactionMechanismGenerator/RMG-Py/tree/main/test/rmgpy). We currently use [pytest](https://docs.pytest.org/en/stable/) for our unit testing.

* If you wrote a significant new feature, please add a regression test:
    * First, please create an input file that includes your new feature. Ensure `saveEdgeSpecies` is set to `True` so that the edge model will also be saved to a file. If applicable, include diverse sets of input conditions to test compatibility with other features.
    * Generate the reaction mechanism corresponding to the input file. Ensure that the simulation does not take more than 15 minutes maximum. You can reduce simulation times in multiple ways, e.g. by increasing the `toleranceMoveToCore` flag.
    * In the `test/regression` folder, create a new folder with a relevant name, and copy the RMG-Py simulation input file in this folder. Include this new folder in your PR.
    * In `.github/workflows/CI.yml`, edit the two lists of regression tests in the `Regression Tests - Execution` and `Regression Tests - Compare to Baseline` steps to add the name of your folder. Be sure to follow BASH syntax.

    * > Warning: This will __fail__ CI because of directory not found errors. This is because the baseline files used for comparison in the regression tests do not exist yet. Your PR will need to be merged by bypassing branch protection restrictions.


### **Do you want to contribute to the documentation?**

* Documentation is [hosted here](http://reactionmechanismgenerator.github.io/RMG-Py/) using [Sphinx](https://www.sphinx-doc.org/en/master/). 

* The live version of the documentation is hosted on the `gh-pages` branch which is updated upon pushes to the `main` branch of RMG-Py.

* To add new documentation, create or modify `.rst` (reStructuredText) files under the `documentation` directory. For a primer on how to write `.rst` markup, please [check out the Sphinx documentation.](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html)

* Please test the documentation on a local build (e.g., via `make html` in the `documentation` directory) before pushing changes.

### **Do you have questions?**

* Email us at rmg_dev@mit.edu.

###  **Best practices for PRs**

* Rebase to the main branch before working, to avoid merge conflicts.

* Keep PRs small and aim to merge quickly.

* Commits should be specific and as small as required; commit messages should be descriptive and as long as required. Commit messages should explain *why* the change is needed. We recommend following [these guidelines.](https://wiki.openstack.org/wiki/GitCommitMessages) 

* Submit a PR only when the code is polished and ready for review. Consider opening a draft PR for work in progress that requires collaborator input.

* Please follow the [PEP8 Python style guide.](https://peps.python.org/pep-0008/)

Thank you!

RMG Developers