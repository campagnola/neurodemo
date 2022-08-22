REM start and run the neurondemo program

SET thispath=%~dp0
echo %thispath:~0,-1%
set "deact=neurodemo_venv\Scripts\deactivate.bat"
set deactcmd="%thispath%%deact%"
echo %deactcmd%
set "act=neurodemo_venv\Scripts\activate.bat"
set actcmd="%thispath%%act%"

call %actcmd%
python demo.py
call %deactcmd%
