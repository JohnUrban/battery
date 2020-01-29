#!/usr/bin/env bash
set -e

echo "
------------------------------------------------------------------------------
                           Looking for cmake make gcc g++ cc c++ ...
------------------------------------------------------------------------------
"

ok=1
for VAR in cmake make gcc g++ cc c++; do
  set +e
  VAL=$(which ${VAR})
  
    if [ $? -eq 0 ]
    then
        echo "  OK $VAR $VAL"
	##${VAR} --version
    else
        echo "  NOT FOUND: $VAR"
        ok=0
    fi
    set -e
done

if [ $ok -ne 1 ]
then
    echo "Please install items NOT FOUND. Cannot continue"
    exit 1
else
    echo "... All OK: 
	Versions not checked. 
	If issues arise, try the recommended versions:
		GNU Make 3.82
		cmake version 3.10.1
		gcc/g++ (GCC) 6.2.0
		cc/c++ (GCC) 4.8.5"
fi
