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
    git pull git@github.com:ReactionMechanismGenerator/RMG-Py.git main

If you are not using ``ssh`` (yet), it is still possible to pull the code using
``https``::
    
    git pull https://github.com/ReactionMechanismGenerator/RMG-Py.git main

We also recommend updating the RMG-database regularly.  The repo itself can be found at https://github.com/ReactionMechanismGenerator/RMG-database/ ::

    cd RMG-database
    git pull git@github.com:ReactionMechanismGenerator/RMG-database.git main

Again, for those (still) using ``https``, the command is instead::

    git pull https://github.com/ReactionMechanismGenerator/RMG-database.git main

We also recommend that the RMS julia package is updated::

    julia -e 'using Pkg; Pkg.update("ReactionMechanismSimulator")'

For more information about how to use the Git workflow to make changes to the source code, please
refer to the handy `Git Tutorial <https://git-scm.com/docs/gittutorial>`_

For information on updating your local repository from ``https`` to ``ssh``,
please see `Managing remote repositories
<https://docs.github.com/en/get-started/git-basics/managing-remote-repositories>`_
