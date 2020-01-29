#!/usr/bin/env bash

ASM_MAP=$1 
READ_MAP=$2
export PATH=~/searchPaths/github/maligner/bin/:${PATH}

QUERY_MISS_PENALTY=3.0
REF_MISS_PENALTY=3.0
QUERY_MAX_MISSES=5
REF_MAX_MISSES=5
SD_RATE=0.05
MIN_SD=750
MAX_SCORE_PER_INNER_CHUNK=1.0
MAX_ALIGNMENTS_PER_QUERY=1


B1=`basename ${ASM_MAP}`
B2=`basename ${READ_MAP}`
OPRE=${B1}_${B2}
maligner_vd \
  -o ${OPRE} \
  -q $QUERY_MISS_PENALTY \
  -r $REF_MISS_PENALTY \
  --query-max-misses $QUERY_MAX_MISSES \
  --ref-max-misses $REF_MAX_MISSES \
  --max-score-per-inner-chunk $MAX_SCORE_PER_INNER_CHUNK \
  --sd-rate $SD_RATE \
  --min-sd $MIN_SD \
  --max-alignments $MAX_ALIGNMENTS_PER_QUERY \
  ${READ_MAP} \
  ${ASM_MAP}

