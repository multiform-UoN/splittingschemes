envname="venv_splitting"
python="python3.8"

rm -rf $envname

$python -m venv $envname

source $envname/bin/activate

$python -m pip install -r requirements.txt

deactivate

clear
echo ""
echo "use the command 'source $envname/bin/activate' in this directory to source the virtual environment!!!"
echo ""
