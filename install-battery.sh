#!/usr/bin/env bash
set -e

function help {
echo "
####################################################################


      ---------------------------------------------------
     |                                     |             |
     |                                     |             |_
(-)  |               BATTERY               |              _| (+)
     |                                     |             |
     |                                     |             |
      ---------------------------------------------------


  Usage: 
	bash install-battery.sh [rfbh]
	-r   Tells this to run the installation. This simply bypasses the catchall help message and exiting. It can be dropped if other options provided.
	-f   To force-install everything anew. This script will not try to re-install a program it already installed on a previous attempt given it finished installation successfully.
	-b   To ignore Boost PreReq fails.
	-h   Help message.
  
  If installation of a program fails, this script will stop. It will not try to install the remaining problems.
  After you figure out the source of the issue, it is a good practice to \"rm -r\" the directory containing the failed installation attempt, before re-running.
####################################################################

  " ## TODO - have each install program remove previous attempt before trying anew.
}

##############################################################################
## DEFAULTS
##############################################################################
HELP=false
FORCE=false
FORCEBOOST=false
RUN=false #Actually, this variable is not used. I just have "-r" in place to force through zero argument help-catchall below.

##############################################################################
## GET OPTS
##############################################################################
while getopts "rfbh" arg; do
    case $arg in
        r) RUN=true;;
        f) FORCE=true;;
        b) FORCEBOOST=true;;
        h) HELP=true; help; exit;;
        *) help; exit;;
    esac
done

## HELP CATCHALL
if [ $# -eq 0 ]; then help; exit; fi




############## EXECUTE #########################################





echo "
####################################################################


      ---------------------------------------------------
     |                                     |             |
     |                                     |             |_
(-)  |               BATTERY               |              _| (+)
     |                                     |             |
     |                                     |             |
      ---------------------------------------------------


####################################################################
"
if ${FORCE}; then FORCE="force"; fi
if ${FORCEBOOST}; then FORCEBOOST=force; fi

cd third_party
bash check-prerequisites.sh ${FORCEBOOST}
bash install-all.sh ${FORCE}
cd ../
MAIN=${PWD}
TEMPLATE=configs/battery_paths_template.cfg
awk -v "BATTERYPATH=${MAIN}" '{sub(/^MAIN=\/path\/to\/battery/,"MAIN="BATTERYPATH); print}' templates/battery_paths.cfg > configs/battery_paths.cfg

## TODO
## CHECK EVERYTHNG IN PATH




## LAST STEP
ln -s src/run_battery.sh battery






echo "
####################################################################


      ---------------------------------------------------
     |                                     |             |
     |                                     |             |_
(-)  |               BATTERY               |              _| (+)
     |                                     |             |
     |                                     |             |
      ---------------------------------------------------


		INSTALLATION COMPLETE

	Add to PATH:
		${MAIN}


	Usage:
		battery setup.cfg 
		battery setup.cfg [battery_paths.cfg]
		battery setup.cfg ${MAIN}/configs/battery_paths.cfg

	Details:
		setup.cfg 
		Specifies locations of the datasets to be used in the analyses.

		An example setup.cfg file is present in:
			${MAIN}/templates/setup.cfg

		Copy that to your working directory and fill in the locations to the datasets you will be using.


		battery_paths.cfg
		Given a successful installation, it is advisable to use the battery_paths.cfg file created in the following location:
			${MAIN}/configs/battery_paths.cfg

		It is used by default when no argument is given to the soft-linked battery script in:
			${MAIN}/battery

		You can also copy that battery_paths.cfg to whatever location you see fit, and provide the path to it there as the 2nd argument.
		The MAIN variable in batter_paths.cfg should be pre-filled with the path to your battery installation.
		The paths to all other software needed will be specified thereafter.

		If you are using different versions of various software, change the path to it manually in battery_paths.cfg



	Data preparation:
	- Shortread
		- takes in paired-end reads
		- the dataset can be any size
		- for most datasets (> 1 million paired reads), LAP will run faster when given a smaller sample (e.g. 1 million paired reads)
			- The same or highly similar rank orders that the full dataset would yield can be expected with small samples
			- See the LAP paper
			- In tests for battery:
				- For Sciara coprophila (~292Mb), highly similar rank orders were seen down to 100,000 paired reads
				- smaller samples were not tested
	- RNA-seq
		- for now, it assumes paired-end short read
		- takes in paired-end reads
		- the dataset can be any size

	- Long Read
		- Can be ONT, PacBio, or other long read technology 
			- provided in fastq 
		- To use long reads with REAPR, FRCbam, and ALE
			- convert to paired-end 
			- helper scripts are supplied to do so
			- Caution/things to consider:
				- TODO
				- Lessons learned when developing this usage

	- Optical/Electronic/Other recognition site maps
		- BioNano: 
			- BNX files need to be converted to the rmap format
			- Helper scripts provided
			- Need to split the smoothed bionano maps up into smaller batches to run in parallel
			- TODO:
				- Some or all of this can be automated as part of the pipelines
		- NabSys
			- Raw Read or NabTig files need to be converted to the correct rmap format
			- May need to then break that up into smaller batches
			- Helper scripts provided
	- BLAST
		- Already-known high quality sequences 
			- e.g. Sanger sequencing of clones, ESTs, etc from the past few decades
			- sequences you are confident in
			- in the future these might be other assemblies or pieces of other assemblies you are confident in
				- using Minimap2
		- Transcriptome assembly
			- Currently uses BLAST
			- Have found that Minimap2 gives highly consistent results and is massively faster
			- TODO - switch it up or provide the option
		- Transcriptome and peptides of 2 different related species
			- Uses blastn, tblastx, tblastn
				- Results from each are perfectly or nearly perfectly correlated
				- In future, just do:
					- tblastx if only transcripts given
					- tblastn if only peptides given (or if both peptides and transcripts given)
			- Currently, good parameters for Minimap2 not found for transcripts
				- Anything under 50-60% identity is lost
				- A tblastn-like usage for Minimap2 would be incredible

####################################################################
"
