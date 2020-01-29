#!/bin/bash

samtools merge --threads $P combined.bam $PBBAM $ONTBAM
