#LINUX
#wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz

#MACOS
#wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-macosx-intel.tar.gz

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
  tar -xzf hmmer-3.1b2-linux-intel-x86_64.tar.gz
elif [ $OS == "mac" ]; then
  tar -xzf hmmer-3.1b2-macosx-intel.tar.gz
else
  echo "Did not run install process as OS is unknown: ${UNAME}. Please provide one of the following OS as the first and only argument: mac, linux."; 
  exit; 
fi
