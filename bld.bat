:: Compile RMG
mingw32-make -j%CPU_COUNT%

:: Install RMG
%PYTHON% setup.py install

:: lazy "install" of everything in our 'external' folder.
:: most of which should probably be elsewhere
mkdir %SP_DIR%\external
xcopy %SRC_DIR%\external %SP_DIR%\external /E /Y
