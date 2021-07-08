set -x   # Show which command is being run

# Download test data

wget -q http://yt-project.org/data/enzo_tiny_cosmology.tar.gz
tar -zxvf enzo_tiny_cosmology.tar.gz
wget -q http://yt-project.org/data/GasSloshingLowRes.tar.gz
tar -zxvf GasSloshingLowRes.tar.gz
wget -q http://yt-project.org/data/FIRE_M12i_ref11.tar.gz
tar -zxvf FIRE_M12i_ref11.tar.gz

# Download answers

wget http://hea-www.cfa.harvard.edu/~jzuhone/${ANSWER_VER}.tar.gz
tar -zxvf ${ANSWER_VER}.tar.gz		

# Set location of yt test data
mkdir -p $HOME/.config/yt
echo "[yt]" > $HOME/.config/yt/yt.toml
echo "test_data_dir = ${PWD}" >> $HOME/.config/yt/yt.toml
cat $HOME/.config/yt/yt.toml

# Install dependencies using conda

conda install --yes numpy pytest pip h5py astropy sherpa cython scipy
conda install --yes -c jzuhone soxs

