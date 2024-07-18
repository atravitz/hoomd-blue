#!/bin/bash

conda env create --file=environment.yml --name popbd_env
conda activate popbd_env

rm -rf build
mkdir build
cd build/

# cd ../hoomd
# ln -s ../correlator-hoomd-plugin/correlator .
# cd ../build

cmake ../ -DCMAKE_CXX_FLAGS=-march=native -DCMAKE_C_FLAGS=-march=native -DCMAKE_INSTALL_PREFIX=`python3 -c "import site; print(site.getsitepackages()[0])"`

make -j4
make install
