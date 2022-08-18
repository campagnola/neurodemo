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
# python3.9 -m venv $ENVNAME || exit 1
py -m venv $ENVNAME || exit 1
# source $ENVNAME/bin/activate || exit 1
source $ENVNAME/Scripts/Activate
# pip3 install --upgrade pip  # be sure pip is up to date in the new env.py

python -m pip install wheel  # seems to be missing (note singular)
python -m pip install cython
# # if requirements.txt is not present, create:
# # pip install pipreqs
# # pipreqs
#
# #Then:
#
python -m pip install -r requirements_local.txt
source $ENVNAME/Scripts/Activate
python --version
python setup_local.py develop
