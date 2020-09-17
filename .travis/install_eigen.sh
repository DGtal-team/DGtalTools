#!/bin/bash


#
# Local install of Eigen on linux system
#

export EIGEN_ROOT="/usr/local/"

if [ $TRAVIS_OS_NAME == linux ];
then
    cd deps
    wget https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.bz2
    bunzip2 eigen-3.3.7.tar.bz2
    tar xf eigen-3.3.7.tar
    

    cd eigen-3.3.7
    mkdir build ; cd build

    cmake .. -DCMAKE_INSTALL_PREFIX="${SRC_DIR}/deps/local"
    make && make install && cd ${SRC_DIR} && export EIGEN_ROOT="$SRC_DIR/deps/local"

fi
