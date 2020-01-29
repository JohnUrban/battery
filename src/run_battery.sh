#!/bin/bash
set -e
echo $0


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
	battery setup.cfg [battery_paths.cfg]

	battery c:[p:aslmrbBh]
	-c Required: Provide path to setup.cfg
	-p Optional: Provide path to battery_paths.cfg. If not provided, default is used.
	-a Optional: Peform all analyses: short, long, maligner, rnaseq, blast, busco. This is same as -s -l -m -r -b -B.
	-s Optional: Perform shortread analysis.
	-l Optional: Perform longread analysis.
	-m Optional: Perform maligner analysis of rmaps.
	-r Optional: Perform RNAseq analysis.
	-b Optional: Perform BLAST analyses.
	-B Optional: Perform BUSCOv3 analyses.
	-h Get help message.
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

	- Long Read
		- Can be ONT, PacBio, or other long read technology 
			- provided in fastq 
		- To use long reads with REAPR, FRCbam, and ALE
			- convert to paired-end 
			- helper scripts are supplied to do so
			- Caution/things to consider:
				- TODO
				- Lessons learned when developing this usage

	- RNA-seq
		- for now, it assumes paired-end short read
		- takes in paired-end reads
		- the dataset can be any size


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
}

if [ $# -eq 0 ]; then help; exit; fi




##############################################################################
## DEFAULTS
##############################################################################

#SETUP=${1}
#BATTERYMAIN=$(dirname ${0})
#if [ $# -gt 1 ]; then BATTERY_PATHS=${2} ;
#else BATTERY_PATHS=${BATTERYMAIN}/configs/battery_paths.cfg ;
#fi

BATTERYMAIN=$(dirname ${0})
BATTERY_PATHS=${BATTERYMAIN}/configs/battery_paths.cfg
DO_SHORT=false
DO_LONG=false
DO_MALIGN=false
DO_RNASEQ=false
DO_BLAST=false
DO_BUSCO=false
HELP=false

##############################################################################
## GET OPTS
##############################################################################
while getopts "c:p:aslmrbBh" arg; do
    case $arg in
	c) SETUP=${OPTARG};;
	p) BATTERY_PATHS=${OPTARG};;
        a) DO_SHORT=true; DO_LONG=true; DO_MALIGN=true; DO_RNASEQ=true; DO_BLAST=true; DO_BUSCO=true;;
        s) DO_SHORT=true;;
        l) DO_LONG=true;;
        m) DO_MALIGN=true;;
        r) DO_RNASEQ=true;;
        b) DO_BLAST=true;;
        B) DO_BUSCO=true;;
        h) HELP=true;;
        *) help; exit;;
    esac
done

##############################################################################
## TRIGGER HELP IF OPTED FOR
##############################################################################
if ${HELP} || [ -z ${SETUP} ] ; then help; exit; fi



##############################################################################
## SET UP ENVIRONMENT GIVEN BATTERY_PATHS AND SETUPCFG
##############################################################################
## NEED TO SOURCE BATTERY PATHS FIRST
source ${BATTERY_PATHS} 

## SETUP SECOND
source ${SETUP} 



##############################################
##############################################
##############################################
################ EXECUTE #####################
NEWDIR=eval_asms
if $RENAME || [ $MINLEN -gt 0 ]; then mkdir -p ${NEWDIR} ; fi

## IF OPTED; rename and/or set min contig length
if $RENAME; then echo renaming....;
 if [ $MINLEN -gt 0 ]; then echo ...also setting min contig length to $MINLEN ; fi
 while read fasta; do
  b=`basename $fasta`
  fasta_name_changer.py -f $fasta -r contig -n --key key-${b}.txt | extractFastxEntries.py --fa --stdin --minlen $MINLEN > ${NEWDIR}/${b}
 done < $ASMFOFN
elif [ $MINLEN -gt 0 ] && [ $RENAME != "true" ]; then echo ...setting min contig length to $MINLEN ... ;
 while read fasta; do
  b=`basename $fasta`
  extractFastxEntries.py --fa -f $fasta --minlen $MINLEN > ${NEWDIR}/${b}
 done < $ASMFOFN
fi

## CREATE AND USE UPDATED FOFN
if $RENAME || [ $MINLEN -gt 0 ]; then
  for f in ${NEWDIR}/*; do
    readlink -f $f;
  done > renamed.fofn
  ASMFOFN=`readlink -f renamed.fofn`
fi

#for f in ${NEWDIR}/*; do
#  readlink -f $f; 
#done > renamed.fofn
#ASMFOFN=`readlink -f renamed.fofn`


## BEGIN LAUNCHING JOBS
if ${DO_SHORT}; then
  echo shortread
  mkdir -p shortread
  cd shortread
  #echo "$SHORTAUTO $ASMFOFN $LR1 $LR2 $R1 $R2 $SHORTCONFIG $SHORTSCRIPTS"
  bash $SHORTAUTO $ASMFOFN $LR1 $LR2 $R1 $R2 $SHORTCONFIG $SHORTSCRIPTS
  cd ../
fi

if ${DO_MALIGN}; then
  echo bionano
  mkdir -p bionano
  cd bionano
  $BIONANORUN $BIONANOCLEAN $BIONANOCONFIG $ASMFOFN $MAPSFOFN $REC_ENZ $REC_SEQ $BIONANOSCRIPTS
  cd ../
fi

if ${DO_LONG}; then
  echo longread
  mkdir -p longread
  cd longread 
  bash $AUTOLR $LRCLEAN $LRCONFIG $ASMFOFN $LRSCRIPTS $ONT $PACBIO $ONT1 $ONT2 $PACBIO1 $PACBIO2
  cd ../
fi


if ${DO_BLAST}; then
  echo blast_analyses
  mkdir -p blast_analyses
  cd blast_analyses

  echo transcriptome
  mkdir -p transcriptome
  cd transcriptome
  $TRANSRUN $TRANSSCRIPTS $TRANSCONFIG $TRANSCLEAN $ASMFOFN $TRANS1 $TRANSNJOBS $TBLASTX $TRANJOBPRE1
  cd ../

  echo dmel
  mkdir -p dmel
  cd dmel
  $TRANSRUN $TRANSSCRIPTS $OTHERSPP_TRANSCONFIG $TRANSCLEAN $ASMFOFN $TRANS2 $TRANSNJOBS $TBLASTX $TRANJOBPRE2
  cd ../

  echo anopheles
  mkdir -p anopheles
  cd anopheles
  $TRANSRUN $TRANSSCRIPTS $OTHERSPP_TRANSCONFIG $TRANSCLEAN $ASMFOFN $TRANS3 $TRANSNJOBS $TBLASTX $TRANJOBPRE3
  cd ../

  echo dmel_peptides
  mkdir -p dmel_peptides
  cd dmel_peptides
  $PEPRUN $PEPSCRIPTS $PEPCONFIG $PEPCLEAN $ASMFOFN $PEP2 $PEPNJOBS $PEPJOBPRE2
  cd ../

  echo anopheles_peptides
  mkdir -p anopheles_peptides
  cd anopheles_peptides
  $PEPRUN $PEPSCRIPTS $PEPCONFIG $PEPCLEAN $ASMFOFN $PEP3 $PEPNJOBS $PEPJOBPRE3
  cd ../

  echo knownseqs
  mkdir -p knownseqs
  cd knownseqs
  $TRANSRUN $TRANSSCRIPTS $KNOWNCONFIG $KNOWNCLEAN $ASMFOFN $KNOWNSEQS $KNOWNNJOBS $KNOWNTBLASTX $KNOWNJOBPRE
  cd ../

  #leave blast_analyses
  cd ../
fi

if ${DO_RNASEQ}; then
  echo rnaseq
  mkdir -p rnaseq
  cd rnaseq
  $RNARUN $RNASCRIPTS $RNACONFIG $RNACLEAN $ASMFOFN $RNAFOFN
  cd ../
fi

if ${DO_BUSCO}; then
  echo buscov3
  mkdir -p buscov3
  cd buscov3
  $BUSCOV3RUN $BUSCOV3SCRIPTS $BUSCOV3CONFIG $BUSCOV3CLEAN $ASMFOFN 
  cd ../
fi
################ EXECUTE #####################
##############################################
##############################################
##############################################





