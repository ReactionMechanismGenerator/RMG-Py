export CC=${PREFIX}/bin/gcc
export CXX=${PREFIX}/bin/g++
export F77=${PREFIX}/bin/gfortran
export F90=${PREFIX}/bin/gfortran
make
$PYTHON setup.py install
