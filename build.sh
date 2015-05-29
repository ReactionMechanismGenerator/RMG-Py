export CC=${PREFIX}/bin/gcc
export CXX=${PREFIX}/bin/g++
export F77=${PREFIX}/bin/gfortran
export F90=${PREFIX}/bin/gfortran
make -j ${CPU_COUNT}
$PYTHON setup.py install

cd ${PREFIX}/lib/python${PY_VER}/site-packages/
mkdir external
cp ${SRC_DIR}/external/* external/ -r 