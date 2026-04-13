set -x

export ATOMDB=$HOME/atomdb
micromamba activate test-env
python -m pytest -vv pyxsim/tests
