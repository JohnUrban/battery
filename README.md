# battery
An automated pipeline for genome assembly evaluations using SLURM.
John Urban (2013-2020)

```
####################################################################


      ---------------------------------------------------
     |                                     |             |
     |                                     |             |_
(-)  |               BATTERY               |              _| (+)
     |                                     |             |
     |                                     |             |
      ---------------------------------------------------


####################################################################
```


Dependencies (not included):
- cmake
- make
- boost
- gcc
- g++
- cc
- c++
- python2.7
- python3
- R
- perl

Dependencies (included):
- sciaraTools
- bamtools (version included: 2.4.1)
- busco 
- augustus
- HMMer (linux or mac)
- ncbi-blast-2.2.30+ (linux or mac)
- bowtie2
- bwa
- minimap2
- blast
- sniffles
- lap
- ale
- reapr
- frcbam
- pilon
- kentTools (linux or mac)
- hisat2
- 


- Python2.7 and Python3
	- Since different systems use the python variable for different versions,
		Battery requires variables in the environment called python2.7 and python3 
		to allow python scripts corresponding to each to work.
	- The headers of some python scripts needed to be changed to reflect this
		- Thus, if your version of battery is not working, it is recommended that you install the versions of those programs packaged with Battery
	- Changed all py scripts in buscoV3 to have python3 in their headers (some already did)
	- Lap -- calc_prob.py python2.7

		
- Python libraries
	- biopython - only for parsing fasta/fastq 
	- numpy
	- argparse
	- collections
	- 


- R

- Perl (majorly for Reapr) - 5.24.1 worked for me
	- File::Basename
	- File::Copy
	- File::Spec
	- File::Spec::Link
	- Getopt::Long
	- List::Util


Recommendations:
- Battery may work best when installing and using the software packages that come with it
- For example, small tweaks were made to some of the software Battery uses
	- ALE
		- A extremely rare issue occurred using ALE with a read that contained a base quality higher than a ceiling ALE anticipates causing ALE to crash
		- The version here allows a slightly higher ceiling to deal with that possibility
		- However, a helper script called enforceQualCeiling.py contained in Battery can also be used to fix an offending fastq file
		- s.t. that base qualities above a given ceiling are changed to be the ceiling
	- LAP
		- There are certain lines in calc_probs.py that were changed and annotated
		- These lines are almost never used during typical uses of LAP, which may be why the bugs they trigger were allowed to exist
	- Maligner
		- By default this program set a pretty low celing on how long the maps can be that you are trying to align
		- The version with Battery has a ceiling that will not likely ever be hit allowing very long maps to align
		- Python scripts have python2.7 in headers instead of python
	- Reapr
		- Reapr is normally packaged with its own cmake and bamtools (the latter depending on the former)
		- For various reasons, it has been more straight forward to remove both from Reapr
		- Instead, the version of Reapr contained in Battery relies on a BamTools that is already installed
		- Since Battery comes with its own version of BamTools, using install-all.sh to install all software packaged with Battery
			- tells Reapr to use the BamTools packaged with Battery
	- BamTools redundancy
	        - BamTools is provided with Battery as a dependency for Augustus
		- However, many of the third_party tools packaged with Battery come with their own BamTools
			- An example being Reapr, discussed above
		- In the future, some of the other tools that carry their own BamTools may also be changed to use the BamTools packaged with Battery instead of their own
		- These include: 
		        - FRC 
		        - Sniffles
		        - BEDtools
		        - FreeBayes (although this program is not yet a part of Battery)











# installation

```
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
        -r   Tells this to run the installation. This simply bypasses the catchall help message and exiting. It can be dropped if other options provid
ed.
        -f   To force-install everything anew. This script will not try to re-install a program it already installed on a previous attempt given it fi
nished installation successfully.
        -b   To ignore Boost PreReq fails.
        -h   Help message.
  
  If installation of a program fails, this script will stop. It will not try to install the remaining problems.
  After you figure out the source of the issue, it is a good practice to \"rm -r\" the directory containing the failed installation attempt, before re
-running.
####################################################################
```




# usage

```
####################################################################


      ---------------------------------------------------
     |                                     |             |
     |                                     |             |_
(-)  |               BATTERY               |              _| (+)
     |                                     |             |
     |                                     |             |
      ---------------------------------------------------


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


```


If you find these tools useful for your own work, please cite the fungus fly genome preprint:
-------------------------------------------
Urban JM, Foulk MS, Bliss JE, Coleman CM, Lu N, Mazloom R, Brown SJ, Spradling AC, Gerbi SA. 2020. 

Single-molecule sequencing of long DNA molecules allows high contiguity de novo genome assembly for the fungus fly, Sciara coprophila. 

bioRxiv 2020.02.24.963009.

https://www.biorxiv.org/content/10.1101/2020.02.24.963009v1
