#!/bin/bash

## One can obtain a merged bnx file by:
## RefAligner -bnx -merge -i /Path/to/RawMolecules1.bnx -i /Path/to/RawMolecules2.bnx -i /Path/to/RawMolecules_N.bnx -o ${PRE} -minSNR 2.75 -minlen 150 -minsites 8 -MaxIntensity 0.6

##############################################################################
## FUNCTION: HELP
##############################################################################

function help {
    echo "
        Usage: ${0} -B:N:o:d:n:h
        -B with argument = input BNG maps (BNX for bionano) -- This arg is mutually exclusive with -N.
        -N with argument = input NabSys maps (nabtigs/nabreads for NabSys). This arg is mutually exclusive with -B.
        -o with argument = output prefix (default: maps)
        -d with argument = name of directory to store split up map locations in (default: splitmaps (subdir in working dir))
        -n with argument = Number of maps per file when splitting up the single large file. (default: 25000).
        -h help - returns this message; also returns this when no arguments given

	The split command used herein only works on Linux distributions in my experience. (-d not supported on MacOS)
	Therefore, if this is used on MacOS - expect an error message about split.
"
}

##############################################################################
## TRIGGER HELP IF NO ARGS
##############################################################################
if [ $# -eq 0 ]; then help; exit; fi



##############################################################################
## DEFAULTS
##############################################################################
## DEFAULTS
PRE=maps
SPLITDIR=splitmaps
NMAPSPERFILE=25000
HELP=false

##############################################################################
## GET OPTS
##############################################################################
while getopts "B:N:o:d:n:h" arg; do
    case $arg in
        B) BNX=$OPTARG;;
        N) NAB=$OPTARG;;
        o) PRE=$OPTARG;;
        d) SPLITDIR=$OPTARG;;
        n) NMAPSPERFILE=$OPTARG;;
        h) HELP=true;;
        *) help; exit;;
    esac
done


##############################################################################
## TRIGGER HELP IF OPTED FOR
##############################################################################
if ${HELP}; then help; exit; fi
if [ ! -z $BNX ] && [ ! -z $NAB ]; then echo "-B and -N are mutually exclusive."; help; exit; fi


if [ ! -z $BNX ]; then bnx2maligner.py -b ${BNX} > ${PRE}.maps; 
elif [ ! -z $NAB ]; then nabsys2maligner.py -f ${NAB} > ${PRE}.maps; 
fi

smooth_maps_file -m 1000 ${PRE}.maps > ${PRE}.smoothed.maps
mkdir -p ${SPLITDIR} && cd ${SPLITDIR}
split -l ${NMAPSPERFILE} -d ../${PRE}.smoothed.maps ${PRE}.smoothed.maps
