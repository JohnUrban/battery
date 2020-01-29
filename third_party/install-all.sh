#!/bin/bash
set -e


## gcc and g++ must be present. Recommend version 6.2. Problems may occur for versions < 5.
## OS information is obtained by the specific install-*sh scripts that need know to know that.

## Boost must be present - for several installs to work. Recommend 1.55.0. 
	# Pertinent variables to check in your environment are:
	#	CPATH LD_LIBRARY_PATH LIBRARY_PATH LD_RUN_PATH BOOST BOOST_DIR BOOST_ROOT BOOST_INCLUDEDIR
	# Try running "check-for-boost.sh" in this directory for more information.
	# Run the install-boost.sh script if need be
	# Will need to make or prepend all relevant variables with location of boost if they are not available already
	# Systems w/ module load boost/1.55.0 Likely do this for you -- so make sure to run that before installation if available (try check-for-boost.sh after that)

function run { 
export PROG=$(echo $PROG + 1 | bc)
if [ -f .INSTALLED_${1} ]; then
  echo "

------------------------------------------------------------------------------

                DETECTED INSTALLATION OF ${PROG}/${NPROGS} : ${1} 
                MOVING ON......

------------------------------------------------------------------------------

"
else
  echo "

------------------------------------------------------------------------------

                INSTALLING ${PROG}/${NPROGS} : ${1} 
			${2} ${3}

------------------------------------------------------------------------------

"
  bash ${2} ${3}
  touch .INSTALLED_${1}
  echo "

------------------------------------------------------------------------------

                FINISHED INSTALLING ${PROG}/${NPROGS} : ${1} 
			${2} ${3}

------------------------------------------------------------------------------

"
fi
}

function clean {
  echo "

------------------------------------------------------------------------------

                FORCE INSTALL OPTED FOR : CLEANING OUT PREVIOUS ATTEMPT

------------------------------------------------------------------------------

"
 rm .INSTALLED_*
 for DIR in */; do
   if [ -d ${DIR} ]; then echo "Removing ${DIR} ....."; rm -rf ${DIR}; fi
 done
}

function help {
  echo "Run: bash install-all.sh

  This script will not try to re-install a program it already installed on a previous attempt given it finished installation successfully.
  To force it to install everything anew, do: 
	bash install-all.sh force
  
  If installation of a program fails, this script will stop. It will not try to install the remaining problems.
  After you figure out the source of the issue, it is a good practice to \"rm -r\" the directory containing the failed installation attempt, before re-running.
  " ## TODO - have each install program remove previous attempt before trying anew.
}


NPROGS=21
PROG=0

## OPTIONAL: CLEAN AND INSTALL
if [ $# -eq 1 ]; then if [ $1 == "force" ]; then clean; fi; fi
if [ $# -eq 1 ]; then if [ $1 == "help" ] || [ $1 == "h" ] || [ $1 == "-help" ] || [ $1 == "--help" ] || [ $1 == "-h" ]; then help; exit; fi; fi


## INSTALLS
run BamTools install-bamtools.sh
run Augustus install-augustus.sh ./bamtools  	## Depends on BamTools install
run REAPR install-reapr.sh ./bamtools		## Depends on BamTools install

## The following install independently
run ALE install-ale.sh
run BEDtools install-bedtools.sh
run BLAST install-blast.sh
run Bowtie2 install-bowtie2.sh
run BUSCOv1 install-buscov1.sh
run BUSCOv3 install-buscov3.sh
run BWA install-bwa.sh
run FRCbam install-frc.sh
run HiSat2 install-hisat2.sh
run HMMer install-hmmer.sh
run kentTools install-kentTools.sh
run LAP install-lap.sh
run Maligner install-maligner.sh
run Minimap2 install-minimap2.sh
run PicardTools install-picard.sh
run Pilon install-pilon.sh
run SAMtools install-samtools.sh
run Sniffles install-sniffles.sh
