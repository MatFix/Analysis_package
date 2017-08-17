SETLOCAL
set LOGFILE_DATE=%DATE:~6,4%.%DATE:~3,2%.%DATE:~0,2%
set LOGFILE_TIME=%TIME:~0,2%.%TIME:~3,2%
set LOGFILE=log-%LOGFILE_DATE%-%LOGFILE_TIME%.log

>%LOGFILE% 2>&1 (
	@SET shutDown="No"
	ECHO "Running file:"
	ECHO.
	"%CD%\nanodot_thermal_Gaussian1.mx3"
	ECHO.
	timeout /t 1500
	ECHO "Running file:"
	ECHO.
	"%CD%\nanodot_thermal_Gaussian2.mx3"
	ECHO. 
	ECHO Run complete.
	TIMEOUT /T 600
	ECHO.
)
	
if %shutDown%=="Yes" (
	ECHO Turning off the computer, press CTRL+C to abort
	TIMEOUT /T 60
	shutdown -s -f -t 0
)