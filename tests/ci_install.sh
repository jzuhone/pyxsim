set -x   # Show which command is being run

if [[ ${mode} == "testing" ]]; then

    # Download test data
    
    wget -q http://yt-project.org/data/enzo_tiny_cosmology.tar.gz
    tar -zxf enzo_tiny_cosmology.tar.gz
    wget -q http://yt-project.org/data/GasSloshingLowRes.tar.gz
    tar -zxf GasSloshingLowRes.tar.gz
    wget -q http://yt-project.org/data/FIRE_M12i_ref11.tar.gz
    tar -zxf FIRE_M12i_ref11.tar.gz
    
    # Download answers
    
    wget -q http://hea-www.cfa.harvard.edu/~jzuhone/${ANSWER_VER}.tar.gz
    tar -zxf ${ANSWER_VER}.tar.gz		
    
    # Set location of yt test data
    mkdir -p $HOME/.config/yt
    echo "[yt]" > $HOME/.config/yt/yt.toml
    echo "test_data_dir = \"${PWD}\"" >> $HOME/.config/yt/yt.toml
    cat $HOME/.config/yt/yt.toml
fi 

# Install dependencies using conda

PYVER=`python --version`

conda config --add channels https://conda.anaconda.org/sherpa
conda install --yes numpy pytest pip h5py astropy nose cython scipy yt
pip install soxs==4.0b6

if [[ ${mode} == "testing" ]]; then
  conda install --yes -c sherpa sherpa
fi 

if [[ ${mode} == "wheels" ]]; then
  conda install --yes wheel setuptools
fi 

if [[ ${mode} == "testing" ]]; then
  # Install pyxsim
  python -m pip install -e .
fi