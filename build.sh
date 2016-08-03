export CC=${PREFIX}/bin/gcc
export CXX=${PREFIX}/bin/g++
export F77=${PREFIX}/bin/gfortran
export F90=${PREFIX}/bin/gfortran
make -j${CPU_COUNT}
make QM
$PYTHON setup.py install

# Save version number stored in rmgpy/__init__.py file
$PYTHON -c 'from rmgpy import __version__; print __version__' > ${SRC_DIR}/__conda_version__.txt

# lazy "install" of everything in our 'external' folder.
# most of which should probably be elsewhere
cp -R ${SRC_DIR}/external ${SP_DIR}
