## Needs boost in ENV

## Forcing it to use gcc/g++ instead of cc/c++ b/c at least on this system, easiest way to get augustus and bamtools to use same version
## I've had luck with gcc/g++ 6.2.0

## Note: FRC comes with its own copy of bamtools - so I could just have other tools link to that or frc link to the standard copy of bamtools
##	For now, just leaving as is...
#git clone https://github.com/vezzi/FRC_align.git

tar -xzf FRC_align.tar.gz 
cd FRC_align
mkdir build
cd build
cmake -D CMAKE_C_COMPILER=$(which gcc) -D CMAKE_CXX_COMPILER=$(which g++) -DCMAKE_INSTALL_PREFIX=../ ..
make


