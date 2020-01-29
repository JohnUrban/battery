#rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ./

MAC=kentTools-macos.tar.gz
LINUX=kentTools-linux.tar.gz


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
  tar -xzf ${LINUX}
elif [ $OS == "mac" ]; then
  tar -xzf ${MAC}
else
  echo "Did not run install process as OS is unknown: ${UNAME}. Please provide one of the following OS as the first and only argument: mac, linux."; 
  exit; 
fi
