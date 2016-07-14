.. _updatingSourceCode:

*******************************
Updating the RMG-Py Source Code
*******************************

It is recommended to keep yourself up to date with the latest patches and bug fixes by RMG developers,
which is maintained on the official repository at https://github.com/ReactionMechanismGenerator/RMG-Py/ 
You can view the latest changes by viewing the commits tab on the repository.  
To update your source code, you can "pull" the latest changes from the official repo by typing the following command in the
Command Prompt ::

    cd RMG-Py
    git pull https://github.com/ReactionMechanismGenerator/RMG-Py.git master

We also recommend updating the RMG-database regularly.  The repo itself can be found at https://github.com/ReactionMechanismGenerator/RMG-database/ ::

    cd RMG-database
    git pull https://github.com/ReactionMechanismGenerator/RMG-database.git master

It is also important to update the conda environment regularly to ensure that all the necessary packages are
found in the active environment. ::

    cd RMG-Py
    conda env update -f environment_[operatingSystemName].yml
    activate rmg_env      OR     source activate rmg_env

For more information about how to use the Git workflow to make changes to the source code, please
refer to the handy `Git Tutorial <http://git-scm.com/docs/gittutorial>`_