#!/bin/bash

mkdir ~/doxygen && cd ~/doxygen
wget https://sourceforge.net/projects/doxygen/files/rel-1.8.14/doxygen-1.8.14.src.tar.gz/download -O doxygen-1.8.14.src.tar.gz \
  && tar xzf doxygen-1.8.14.src.tar.gz
cd doxygen-1.8.14 ; mkdir build ; cd build
cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$TRAVIS_BUILD_DIR/doxygen/
ninja install
