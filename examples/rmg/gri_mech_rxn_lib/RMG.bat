@echo off

REM RMG-Py Windows batch script for RMG execution
REM Put me in the directory containing the input file and double-click to run RMG.
REM This assumes that the condition file is called input.py.
REM Output from RMG will be logged to the file RMG.log.

echo Running RMG...
python "..\..\..\rmg.py" input.py
echo RMG job completed.

:end
pause
