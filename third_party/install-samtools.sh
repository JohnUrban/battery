#wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2
#wget https://github.com/samtools/htslib/releases/download/1.3/htslib-1.3.tar.bz2

tar -xzf samtools.tar.gz
cd samtools

MAIN=${PWD} ##samtools top dir

mkdir -p ${MAIN}/localinstall

tar -xjf htslib-1.3.tar.bz2
cd htslib-1.3    
make
make prefix=${MAIN}/localinstall install
cd ../

tar -xjf samtools-1.3.tar.bz2
cd samtools-1.3    
make
make prefix=${MAIN}/localinstall install
cd ../



