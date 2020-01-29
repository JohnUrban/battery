#!/usr/bin/env bash
set -e

## NOTE: Sections of this code are borrowed from REAPR's install.sh script. (Thanks)
echo "
------------------------------------------------------------------------------
                           Checking Perl modules...
------------------------------------------------------------------------------
"

modules_ok=1

for module in File::Basename File::Copy File::Spec File::Spec::Link Getopt::Long List::Util
do
    set +e
    perl -M$module -e 1 2>/dev/null

    if [ $? -eq 0 ]
    then
        echo "  OK $module"
    else
        echo "  NOT FOUND: $module"
        modules_ok=0
    fi
    set -e
done

if [ $modules_ok -ne 1 ]
then
    echo "Some Perl modules were not found - please install them. Cannot continue"
    exit 1
else
    echo "... Perl modules all OK"
fi

