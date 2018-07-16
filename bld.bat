:: Install RMG
mingw32-make install

:: lazy "install" of everything in our 'external' folder.
:: most of which should probably be elsewhere
mkdir %SP_DIR%\external
xcopy %SRC_DIR%\external %SP_DIR%\external /E /Y
