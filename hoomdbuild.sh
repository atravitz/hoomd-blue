#!/bin/bash

# Remove old build directory
# rm -rf build

# Make build directory
# mkdir -p build
cd build 
# Define some variables
export CC=$(which gcc)
export CXX=$(which g++)
echo "Using compilers $($CC --version | head -n 1), $($CXX --version | head -n 1)."

# Some cmake setup
# CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=${SOFTWARE_ROOT}/lib/python"

# Compile against correct python
CMAKE_FLAGS="${CMAKE_FLAGS} -DPYTHON_EXECUTABLE=$(which python)"
PYTHON_LIBRARY_PATH=$(python -c "import distutils.sysconfig as sysconfig; import os; print(os.path.join(sysconfig.get_config_var('LIBDIR'), sysconfig.get_config_var('LDLIBRARY')))")
CMAKE_FLAGS="${CMAKE_FLAGS} -DPYTHON_LIBRARY=${PYTHON_LIBRARY_PATH}"  # -DCMAKE_BUILD_TYPE=Debug"
CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBUILD_HPMC=OFF -DBUILD_METAL=OFF -DBUILD_CGCMM=OFF -DENABLE_MPI=OFF"
# Install to the conda packages path
CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX}/lib/python3.8/site-packages"

# link correlator
# cd ../hoomd 
# ln -s ../correlator-hoomd-plugin/correlator .

# cd ../build

cmake ../ ${CMAKE_FLAGS} 
make -j4  
make install 
