#!/bin/bash

ALE ${BAM} ${REF} ${BASE}.ALE.txt >> ${BASE}.ale.err

if $CLEAN; then bash $SCRIPTS/ale.clean.sh ; fi

