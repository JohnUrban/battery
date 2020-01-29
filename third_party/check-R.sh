#!/usr/bin/env bash
set -e

## NOTE: Sections of this code are borrowed from REAPR's install.sh script. (Thanks)

echo "
------------------------------------------------------------------------------
                           Looking for R...
------------------------------------------------------------------------------
"
if type -P R
then
    echo "... found R OK"
else
    echo "Didn't find R. It needs to be installed and in your path. Cannot continue"
fi


