#wget http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30%2B-universal-macosx.tar.gz

### MAC=ncbi-blast-2.2.30+-universal-macosx.tar.gz
### LINUX=ncbi-blast-2.2.30+-x64-linux.tar.gz

MAC=ncbi-blast-2.2.30+MACOS
LINUX=ncbi-blast-2.2.30+LINUX

## OBTAIN OS INFORMATION
UNAME=`uname`
OS=""
if [ $# -eq 0 ]; then 
  if [ $UNAME == "Darwin" ]; then OS=mac; ## Overwrites user arg if given
  elif [ $UNAME == "Linux" ]; then OS=linux; ## Overwrites user arg if given
  fi
else OS=$1
fi

if [ $OS == "linux" ]; then 
  ### tar -xzf ${LINUX}
  ln -s ${LINUX} ncbi-blast-2.2.30+
elif [ $OS == "mac" ]; then
  ### tar -xzf ${MAC}
  ln -s ${MAC} ncbi-blast-2.2.30+
else
  echo "Did not run install process as OS is unknown: ${UNAME}. Please provide one of the following OS as the first and only argument: mac, linux."; 
  exit; 
fi
