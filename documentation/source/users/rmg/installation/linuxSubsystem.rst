.. _linuxSubsystem:

*****************************************************
Installing RMG in the Linux Subsystem on Windows 10
*****************************************************

Requirements
==============

In order to install the Linux Subsystem, you need to be running Windows 10 build 16215 or later. The build can be
determined from the "About" tab in the system settings (see the "OS build" line).

Installing the Linux Subsystem
===================================

1. Follow the instructions provided by Microsoft to install the Linux subsystem, available `here
   <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_. The basic steps include enabling the Windows Linux
   subsystem from a powershell **run as an administrator**, and downloading the latest LTS version of Ubuntu
   *from the Windows store*. We recommend Ubuntu (for which these instructions were made) over the other Linux
   distributions.

2. Once the Linux subsystem is installed, open a web browser in Windows and go to the
   `Anaconda Python Platform Downloads Page <https://www.anaconda.com/download/#linux>`_. Go to the tab for the
   **Linux Installer**, and **right click** on the download icon for Python 3.7 to copy the link location. Open an Ubuntu
   terminal (type in ``Ubuntu`` into the Windows search bar if you are unsure where to find it) and paste the link
   into the terminal immediately after typing the ``wget`` command, so that your terminal looks like the following: ::

    wget https://repo.anaconda.com/archive/Anaconda3-2019.07-Linux-x86_64.sh

   Your exact link will look similar to the one above, but may be a more recent version of the installer. Execute this
   command in the terminal to begin downloading the installer.

3. Once the Anaconda installer has downloaded, execute the following commands in the Ubuntu terminal, changing the name
   of ``Anaconda3-2019.07-Linux-x86_64.sh`` to match the name of the script you downloaded. ::

    bash Anaconda3-2019.07-Linux-x86_64.sh

   Install the anaconda3 folder inside your home directory (it should be the default location when it asks for a location
   to install). **When prompted to append Anaconda to your PATH, select or type Yes**. When prompted, do NOT install
   Microsoft VSCode. If you are interested in this lightweight IDE then you will want to install this into Windows 10
   instead of inside the linux subsystem.

4. Execute the following commands to make sure that all of the required packages in Ubuntu are also installed: ::

    sudo apt install gcc
    sudo apt install g++
    sudo apt install make
    sudo apt-get install libxrender1

5. Follow the instructions for either the binary (:ref:`anacondaUser`) or source installation (:ref:`anacondaDeveloper`)
   for the Linux Operating system. Follow these instructions from the point directly after installing Anaconda.
