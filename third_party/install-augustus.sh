##wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.2.2.tar.gz

## I've had success with: Boost 1.55.0, zlib 1.2.11, 
## gcc/g++ 6.2 --  gcc (GCC) 6.2.0
## 
## NOTE - "-local" simply means it relies on a local version of BAMtools passed to instead of a global BAMtools hard-coded -- see README-JMU inside augustus dir.
## bamtools bin/ needs to be in PATH
## bamtools lib needs to be in LD_LIBRARY_PATH
## bamtools include needs to be in CPATH
## export PATH=${PWD}/bamtools/bin:${PATH}
## export LD_LIBRARY_PATH=${PWD}/bamtools/lib:${LD_LIBRARY_PATH}
## export CPATH=${PWD}/bamtools/include:${CPATH}


if [ $# -eq 0 ]; then echo "Please provide path to bamtools dir (above bin/lib/include) as first and only argument."; exit 1; fi

BAMTOOLS_PATH=`readlink -f ${1}` ## e.g. ./bamtools

##tar -xzf augustus-3.2.2-local.tar.gz
cd augustus-3.2.2-local
make BAMTOOLS=${BAMTOOLS_PATH}
cd config

## ADD FOLLOWING LINE TO BASH PROFILE AND/OR SCRIPTS
export AUGUSTUS_CONFIG_PATH=${PWD}

