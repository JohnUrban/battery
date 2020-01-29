#PRE=Nabtigs_reads_maps.smoothed.rmap_falcon-1kb.maps
PRE=$1
# HEADER = chr1, start1, end1, chr2, start2, end2, sumMscore, strand1, strand2, queryCoords1, queryCoords2, score1, score2, expectedQueryLength
# Expected query length can be used to help infer deletion or insertion on reference.....

paste <(awk 'OFS="\t" {print $2,$10,$11, $1":"$8"-"$9,$19,$3}' ${PRE}.pfx.aln ) <( awk 'OFS="\t" {print $2,$10,$11, $1":"$8"-"$9, $19, $3, $1, $9}' ${PRE}.sfx.aln ) | \
	awk 'OFS="\t" {if($6=="F") $6="+"; if($6=="R") $6="-"; if($12=="F") $12="+"; if($12=="R") $12="-"; print $1,$2,$3,$7,$8,$9,$13,$5+$11,$6,$12,$4,$10, $5, $11, $14}'

