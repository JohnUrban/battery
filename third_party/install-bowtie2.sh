#wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-macos-x86_64.zip         
#wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip
#unzip bowtie2-2.2.9-linux-x86_64.zip
#unzip bowtie2-2.2.9-macos-x86_64.zip


MAC=bowtie2-2.2.9-macos-x86_64.zip
LINUX=bowtie2-2.2.9-linux-x86_64.zip


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
  unzip ${LINUX}
elif [ $OS == "mac" ]; then
  unzip ${MAC}
else
  echo "Did not run install process as OS is unknown: ${UNAME}. Please provide one of the following OS as the first and only argument: mac, linux."; 
  exit; 
fi
