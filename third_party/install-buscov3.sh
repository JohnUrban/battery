

mkdir -p localpy
MAIN=$PWD
PRE=${MAIN}/localpy



tar -xzf buscov3.tar.gz 
cd buscov3/
python setup.py install --prefix="${PRE}"

## next 2 cmds (head and echo) will create and tailor config.ini file for my Oscar acct
head -n 45 config/config.ini.default > config/config.ini

echo "
[tblastn]
# path to tblastn
path = ${MAIN}/ncbi-blast-2.2.30+/bin/

[makeblastdb]
# path to makeblastdb
path = ${MAIN}/ncbi-blast-2.2.30+/bin/

[augustus]
# path to augustus
path = ${MAIN}/augustus-3.2.2-local/bin/

[etraining]
# path to augustus etraining
path = ${MAIN}/augustus-3.2.2-local/bin/

# path to augustus perl scripts, redeclare it for each new script
[gff2gbSmallDNA.pl]
path = ${MAIN}/augustus-3.2.2-local/scripts/ 
[new_species.pl]
path = ${MAIN}/augustus-3.2.2-local/scripts/ 
[optimize_augustus.pl]
path = ${MAIN}/augustus-3.2.2-local/scripts/

[hmmsearch]
# path to HMMsearch executable -- ASSUMES ONLY MACOS OR LINUX VERSION WAS UNTARRED
path = ${MAIN}/hmmer-3.1b2/binaries/

[Rscript]
# path to Rscript, if you wish to use the plot tool
path = $(dirname $(which Rscript))
" >> config/config.ini
