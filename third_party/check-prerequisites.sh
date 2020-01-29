#!/usr/bin/env bash
set -e

DIR=$(dirname ${0})

echo "
------------------------------------------------------------------------------
                           Checking prerequisites
------------------------------------------------------------------------------
"

bash ${DIR}/check-python.sh
bash ${DIR}/check-R.sh
bash ${DIR}/check-perl.sh
bash ${DIR}/check-for-boost.sh ${1}


echo "
------------------------------------------------------------------------------
                           Finished checking prerequisites
------------------------------------------------------------------------------
"
