envname="venv"
python="python3.11"

rm -rf $envname

$python -m venv $envname

source $envname/bin/activate

pip install --upgrade pip

python -m pip install --no-cache -r requirements.txt

deactivate

clear
echo ""
echo ">>> Use the command 'source $envname/bin/activate' in this directory to source the virtual environment!!!"
echo ""
