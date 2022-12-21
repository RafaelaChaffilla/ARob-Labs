set MATLAB=C:\Program Files\Polyspace\R2019a

cd .

if "%1"=="" ("C:\PROGRA~1\POLYSP~1\R2019a\bin\win64\gmake"  -f ARDroneTTSim.mk all) else ("C:\PROGRA~1\POLYSP~1\R2019a\bin\win64\gmake"  -f ARDroneTTSim.mk %1)
@if errorlevel 1 goto error_exit

exit 0

:error_exit
echo The make command returned an error of %errorlevel%
exit 1
