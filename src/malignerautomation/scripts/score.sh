#!/bin/bash

## SUM up scores from ALL alignments, score = 0 if no alignments
tail -n +2 ${ALL} | awk -v "SUM=0" '{SUM += $19} END {print SUM}' > score.txt

