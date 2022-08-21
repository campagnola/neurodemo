set -e # force failure if anyting fails in script - ensuring completion
set -o errexit
ENVNAME="neurodemo_venv"
if [ -d $ENVNAME ]
then
    echo "Removing previous environment: $ENVNAME"
    set +e
    # rsync -aR --remove-source-files $ENVNAME ~/.Trash/ || exit 1
    set -e
    rm -Rf $ENVNAME 
else
    echo "No previous environment to remove."
fi

source $ENVNAME/Scripts/Activate
# pip3 install --upgrade pip  # be sure pip is up to date in the new env.py

python -m pip install wheel  # seems to be missing (note singular)
# python -m pip install cython
#
python -m pip install -U -r requirements_local.txt
# source $ENVNAME/Scripts/Activate
python --version
python setup_local.py develop
