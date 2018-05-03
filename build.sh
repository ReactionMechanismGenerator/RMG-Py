export CC=${GCC}
export CXX=${PREFIX}/bin/g++
export F77=${PREFIX}/bin/gfortran
export F90=${PREFIX}/bin/gfortran
make -j${CPU_COUNT}
make check

$PYTHON ${SRC_DIR}/setup.py install
# lazy "install" of everything in our 'external' folder.
# most of which should probably be elsewhere
cp -R ${SRC_DIR}/external ${SP_DIR}
