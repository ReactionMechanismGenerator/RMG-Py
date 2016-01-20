set CC=gcc
set CXX=g++
set F77=gfortran
set F90=gfortran


mingw32-make -j%CPU_COUNT%
mingw32-make QM

%PYTHON% setup.py install

:: Save version number stored in rmgpy/__init__.py file
%PYTHON% -c "from rmgpy import __version__; print __version__" > %SRC_DIR%\__conda_version__.txt

:: lazy "install" of everything in our 'external' folder.
:: most of which should probably be elsewhere
mkdir %SP_DIR%\external
xcopy %SRC_DIR%\external %SP_DIR%\external /E /Y
