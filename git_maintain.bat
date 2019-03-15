@echo OFF
REM git init

:Start
set /p branch=What branch to access?
git checkout %branch%
echo.
REM git pull
REM echo.

REM echo. Make changes to desired file(s)
REM pause
REM echo.
git status
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
IF /I %cont%==y (GOTO Commit)
set /p cont=Reset HEAD?
IF /I %cont%==y (git reset HEAD)
GOTO Tracking

:Commit
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