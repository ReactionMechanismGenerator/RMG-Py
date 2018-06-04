.. _faq:

******************
FAQ collection
******************


* Got an error of ``Segmentation fault:11`` after installing RMG on my machine?

	**Segmentation fault** is a typical error in C code, caused by a program trying to read or write an illegal memory location, i.e. one it is not allowed to access. The most common cause in RMG is a conflict between two different versions of a shared library.  RMG has some dependencies which are written in C++, e.g. rdkit, openbabel. If you compile one of these with a different version of some compiler library, or you compile RMG using one version and run it with another, you will often get a Segmentation fault. Chances are those packages are not up to date, or maybe your environmental variable ``PATH`` is messed up so that the wrong version of something is being found. Please see one example from a user having same `Segmentation fault issue <https://github.com/ReactionMechanismGenerator/RMG-website/issues/125>`_.

* How can I install RMG-Py without Anaconda?

	Usually we don't recommend installing RMG-Py without Anaconda because it takes longer and is easier to get trouble with package management. But one still can try direct installation on Linux or MacOS by following :ref:`Linux instruction<linux>` or :ref:`MacOS instruction<macos>`. The RMG team does not use this install approach internally any more, so these instructions are not actively maintained.

* Windows binary installation gives ``WindowsError: [Error 5]``?
	
	Error 5 is access is denied, so this is either a permissions error, or an issue with the Windows file lock. `These posts <https://github.com/conda/conda/issues/708>`_ suggest rebooting the computer (in case it's a file lock), and running the anaconda prompt, from which you run ``conda create -c rmg --name rmg_env rmg rmgdatabase``, as an administrator (in case it's a permissions error). Please checkout one example from a user having `Windows binary installation issue <https://github.com/ReactionMechanismGenerator/RMG-Py/issues/779>`_.
	
* I get something like ``IOError: [Errno 13] Permission denied: 'C:\\RMG.log'`` 

	You do not have permission to write to the log file. Try running the RMG from a different folder that you do have write permission to, such as within your user's documents directory, or else try running the command prompt as an Administrator (so that you have write permission everywhere). See for example `issue #817 <https://github.com/ReactionMechanismGenerator/RMG-Py/issues/817>`_.

If you have any other errors please report them by opening an `issue <https://github.com/ReactionMechanismGenerator/RMG-Py/issues?q=is%3Aissue>`_, and for general questions ask in the `RMG-Py chat room <https://gitter.im/ReactionMechanismGenerator/RMG-Py>`_.
