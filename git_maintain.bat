@echo OFF
git init

:Start
set /p branch=What branch to access?
git checkout %branch%
echo.
git pull
echo.

echo. Make changes to desired file(s)
pause
echo.
:Tracking
set track=""
set /p track=What file was changed?
IF %track%=="" (GOTO Committing)
git add %track%
GOTO Tracking

:Committing
git status
echo.
set /p cont=Are the modified files displayed correctly?
IF /I %cont%==n (GOTO Tracking)
set /p msg=Commit message (Use " "):
git commit -m %msg%
echo.

set /p Push=Upload to remote database?
IF /I %Push%==y (git push -u origin %branch%)
echo.
echo.
set /p cont=Still working?
IF /I %cont%==y (GOTO Tracking)

echo Exiting...
TIMEOUT /T 2
exit