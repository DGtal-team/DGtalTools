#!/bin/bash

#
# Note: gmp and boost already installed
#
brew install ccache
export PATH="/usr/local/opt/ccache/libexec:$PATH"

brew update > brew.log
#brew install qt5 doxygen homebrew/science/hdf5 graphviz graphicsmagick fftw eigen
brew install qt5 graphicsmagick fftw eigen
# Explicit install of libqglviewer
# Explicit install of libqglviewer
brew tap DGtal-team/homebrew-dgtal
brew install libqglviewer

## Temporary HDF5 build issue
export BTYPE="$BTYPE -DWITH_HDF5=false" && echo "Disabling HDF5 on MacOS";
export BTYPE="$BTYPE -DWITH_QT5=true -DCMAKE_PREFIX_PATH=$(brew --prefix qt5)" && echo "Forcing Qt5 on MacOS";
