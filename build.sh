# Compile RMG
make -j${CPU_COUNT}

# Install RMG
$PYTHON setup.py install

# lazy "install" of everything in our 'external' folder.
# most of which should probably be elsewhere
cp -R ${SRC_DIR}/external ${SP_DIR}
