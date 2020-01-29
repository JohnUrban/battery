#!/usr/bin/env bash
set -e

echo "
------------------------------------------------------------------------------
                           Checking for Python2.7... and required modules.
------------------------------------------------------------------------------
"

PY2=$(which python2.7)
if [ -z ${PY2} ]; then echo "Could not find python2.7: you may need to softlink python2.7-->python in same directory python is found."; exit 1;
else echo ${PY2}
fi


modules_ok=1

for module in Bio numpy random
do
    set +e
    python2.7 -c "import ${module}" 2>/dev/null

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
    echo "Some Python2.7 modules were not found - please install them. Cannot continue"
    exit 1
else
    echo "... Python2.7 modules all OK"
fi


echo "
------------------------------------------------------------------------------
                           Checking for Python3...
------------------------------------------------------------------------------
"
PY3=$(which python3)
if [ -z ${PY3} ]; then echo "Could not find python3: you may need to softlink python3-->python in same directory python is found."; exit 1;
else echo ${PY3}
fi

modules_ok=1

for module in sys  ## sys comes default - just using as a place-holder for now.
do
    set +e
    python2.7 -c "import ${module}" 2>/dev/null

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
    echo "Some Python3 modules were not found - please install them. Cannot continue"
    exit 1
else
    echo "... Python3 modules all OK"
fi


