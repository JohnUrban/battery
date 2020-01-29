#!/bin/bash

function contain { echo ${!1} | grep -i ${2} | wc -l; }

FORCE=no
if [ $# -gt 0 ]; then if [ $1 == "force" ]; then FORCE=yes; fi; fi


echo "
------------------------------------------------------------------------------
                           Checking for Boost...
------------------------------------------------------------------------------
"

ERROR_MSG=""
ERRORS=0
## Checking for Boost
for VAR in CPATH LD_LIBRARY_PATH LIBRARY_PATH LD_RUN_PATH BOOST BOOST_DIR BOOST_ROOT BOOST_INCLUDEDIR; do
  CONTAIN=`contain ${VAR} boost`
  if [ $CONTAIN -lt 1 ]; then ERRORS=1; ERROR_MSG+="Not all Boost variables set: ${VAR}. \n"; fi
done

## Check for GCC 5.4 or higher...

if [ ${ERRORS} -gt 0 ] && [ ${FORCE} == "no"  ]; then 
	echo "There were errors detected with respect to the environment."; 
	echo "Errors here will lead to failures of some programs to install and run."
	echo -e ${ERROR_MSG}; 
	exit 1
elif [ ${ERRORS} -gt 0 ] && [ ${FORCE} == "yes" ]; then 
	echo "There were errors detected with respect to the environment."; 
	echo "Errors here will lead to failures of some programs to install and run."
	echo "You have chosen to force-try the installation anyway...."
	echo -e ${ERROR_MSG}; 
	exit 0
else echo "Seems like boost may be present."
fi


