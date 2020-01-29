#!/bin/bash

F1=`ls ${ALN}/*aln | head -n 1`
head -n 1 $F1 > all.bionano.smoothed.maps.aln

for f in ${ALN}/*.smoothed.maps*aln; do
    tail -n +2 $f >> all.bionano.smoothed.maps.aln
done
tail -n +2 all.bionano.smoothed.maps.aln | wc -l > num_alignments.txt
