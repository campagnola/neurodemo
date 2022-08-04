set -e # force failure if anyting fails in script - ensuring completion
set -o errexit
ENVNAME="neurodemo_venv"
if [ -d $ENVNAME ]
then
    echo "Removing previous environment: $ENVNAME"
    set +e
    #rsync -aR --remove-source-files $ENVNAME ~/.Trash/ || exit 1
    set -e
    rm -Rf $ENVNAME
else
    echo "No previous environment to remove."
fi
pip3 install virtualenv
python3 -m virtualenv -p /usr/local/bin/python3.10 $ENVNAME
# python3.10 -m venv $ENVNAME || exit 1
source $ENVNAME/bin/activate || exit 1
python3 -m pip install --upgrade pip

pip3 install --upgrade pip  # be sure pip is up to date in the new env.
pip3 install wheel  # seems to be missing (note singular)
pip3 install cython
# # if requirements.txt is not present, create:
# # pip install pipreqs
# # pipreqs
#
# #Then:
#
pip3 install -r requirements_local.txt
source $ENVNAME/bin/activate
python --version
python setup_local.py develop
