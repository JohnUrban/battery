#!/usr/bin/env bash

##Need: ASM_FASTA, BASE, REC_ENZ, REC_SEQ, MIN_FRAG_SIZE

# outputs
ASM_OUT_PFX=${BASE}.${REC_ENZ}

# convert the asm fasta file to the Maligner maps format and smooth the maps file by merging consecutive fragments that are less than 1kb
make_insilico_map -o $ASM_OUT_PFX $ASM_FASTA $REC_SEQ
smooth_maps_file -m $MIN_FRAG_SIZE ${ASM_OUT_PFX}.maps > ${ASM_OUT_PFX}.smoothed.maps

