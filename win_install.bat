REM Install a virtual environment for demo.py
REM Run this under Windows PowerShell
REM
SET thispath=%~dp0
echo %thispath:~0,-1%
set "deact=neurodemo_venv\Scripts\deactivate.bat"
set deactcmd="%thispath%%deact%"
echo %deactcmd%
set "act=neurodemo_venv\Scripts\activate.bat"
set actcmd="%thispath%%act%"

git pull https://github.com/pbmanis/neurodemo
git checkout PyQt6
call %deactcmd%


del /F/Q neurodemo_venv

which python
python -m pip install virtualenv
python -m virtualenv -p python3.10 neurodemo_venv

call %actcmd%
which python
REM python.exe -m pip install --upgrade pip  REM be sure pip is up to date in the new env - do on command line separately

python -m pip install wheel
python -m pip install cython

python -m pip install -U -r requirements_local.txt
echo "Install of all requirements completed"
call %actcmd%
python --version
python setup_local.py develop

