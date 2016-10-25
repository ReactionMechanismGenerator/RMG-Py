.. _faq:

******************
FAQ collection
******************


* Got an error of ``Segmentation fault:11`` after installing RMG on my machine?

	**Segmentation fault** is a typical error in C code. RMG has some dependencies which are written in C++, e.g. rdkit, openbabel. Chances are those packages are not up to date or maybe environmental variable ``PATH`` is messed up, please see one example from a user having same `Segmentation fault issue <https://github.com/ReactionMechanismGenerator/RMG-website/issues/125>`_.


* How can I install RMG-Py without Anaconda?

	Usually we don't recommend installing RMG-Py without Anaconda because it's taking longer and easier to get trouble in package management and RMG team is not maintaining this install approach any more. But one still can try direct installation on Linux or MacOS by following :ref:`Linux instruction<linux>` or :ref:`MacOS instruction<macos>`. 

* Windows binary installation gives ``WindowsError: [Error 5]``?
	
	Error 5 is access is denied, so this is either a permissions error, or an issue with the Windows file lock. `These posts <https://github.com/conda/conda/issues/708>`_ suggest rebooting the computer (in case it's a file lock), and running the anaconda prompt, from which you run ``conda create -c rmg --name rmg_env rmg rmgdatabase``, as an administrator (in case it's a permissions error). Please checkout one example from a user having `Windows binary installation issue <https://github.com/ReactionMechanismGenerator/RMG-Py/issues/779>`_.
	
