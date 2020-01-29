#!/bin/bash

##git clone git://github.com/pezmaster31/bamtools.git
## 

## Forcing it to use gcc/g++ instead of cc/c++ b/c at least on this system, easiest way to get augustus and bamtools to use same version
## I've had luck with gcc/g++ 6.2.0

tar -xzf bamtools.tar.gz
cd bamtools
mkdir build
cd build
cmake -D CMAKE_C_COMPILER=$(which gcc) -D CMAKE_CXX_COMPILER=$(which g++) -DCMAKE_INSTALL_PREFIX=../ ..
make
cd ..

#mkdir -p ../local/lib ../local/include
#cp -r include ../local/include/bamtools
#cp -r lib ../local/lib/bamtools

