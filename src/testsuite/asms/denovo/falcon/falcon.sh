#!/bin/bash


mkdir $DIR
cp fc_run.cfg $DIR/fc_run.cfg
cd $DIR
fasta2falcon.py -f $READS > reads.fasta
readlink -f reads.fasta > input.fofn



export PATH=~/software/falcon73/FALCON-integrate/fc_env/bin:/usr/local/bin:/gpfs/runtime/opt/python/2.7.3/bin:/bin:/usr/bin #:$PATH
export PYTHONPATH=/users/jurban/software/falcon73/FALCON-integrate/FALCON/:~/software/falcon73/FALCON-integrate/pypeFLOW/:~/software/falcon73/FALCON-integrate/fc_env/lib/python2.7/site-packages/:/users/jurban/software/localpy/lib/python2.7/site-packages
fc_run.py fc_run.cfg /users/jurban/software/falcon73/logging.ini
