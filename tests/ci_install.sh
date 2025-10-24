set -x   # Show which command is being run

if [[ ${mode} == "testing" ]]; then

    # Download test data

    curl -OL http://yt-project.org/data/enzo_tiny_cosmology.tar.gz
    tar -zxf enzo_tiny_cosmology.tar.gz
    curl -OL http://yt-project.org/data/GasSloshingLowRes.tar.gz
    tar -zxf GasSloshingLowRes.tar.gz
    curl -OL http://yt-project.org/data/FIRE_M12i_ref11.tar.gz
    tar -zxf FIRE_M12i_ref11.tar.gz

    # Download answers

    curl -OL http://hea-www.cfa.harvard.edu/~jzuhone/${ANSWER_VER}.tar.gz
    tar -zxf ${ANSWER_VER}.tar.gz

    # Set location of yt test data
    mkdir -p $HOME/.config/yt
    echo "[yt]" > $HOME/.config/yt/yt.toml
    echo "test_data_dir = \"${PWD}\"" >> $HOME/.config/yt/yt.toml
    cat $HOME/.config/yt/yt.toml

    # Set location of SOXS data
    mkdir -p $HOME/.config/soxs
    echo "[soxs]" > $HOME/.config/soxs/soxs.cfg
    echo "soxs_data_dir = ${GITHUB_WORKSPACE}/soxs_data" >> $HOME/.config/soxs/soxs.cfg
    echo "soxs_answer_dir = ${GITHUB_WORKSPACE}/${ANSWER_VER}" >> $HOME/.config/soxs/soxs.cfg
    cat $HOME/.config/soxs/soxs.cfg

fi

# Install dependencies using conda

# Install dependencies using mamba and pip

eval "$(micromamba shell hook --shell bash)"
micromamba shell init --shell bash --root-prefix=~/micromamba
micromamba activate test-env

micromamba install --yes -c conda-forge numpy pytest pip astropy scipy cython h5py yt

git clone https://github.com/lynx-x-ray-observatory/soxs
cd soxs
pip install .
cd ..

if [[ ${mode} == "wheels" ]]; then
  conda install --yes wheel setuptools
fi

if [[ ${mode} == "testing" ]]; then
  # Install pyxsim
  python -m pip install -e .
fi
