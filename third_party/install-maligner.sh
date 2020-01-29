#git clone https://github.com/JohnUrban/maligner.git

tar -xzf maligner-JMU.tar.gz

cd maligner
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -D CMAKE_C_COMPILER=$(which gcc) -D CMAKE_CXX_COMPILER=$(which g++) -DCMAKE_INSTALL_PREFIX=../ .. .. 2>&1 | tee cmake.out
make
make install


## does this need gcc 4.9.2
##cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/gpfs/runtime/opt/gcc/4.9.2/bin/c++ .. 2>&1 | tee cmake.out
