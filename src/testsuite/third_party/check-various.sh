#!/bin/bash
set -e

echo "Checking for fastqSimulate from Canu scripts..."
echo "Found: " $(which fastqSimulate)

echo "Checking for faSize from kentTools..."
echo "Found: " $(which faSize)

echo "Checking for bedtools...."
echo "Found: " $(which bedtools)

echo "Checking for make_insilico_map"
echo "Found: " $(which make_insilico_map)

echo "Checking for smooth_maps_file"
echo "Found: " $(which smooth_maps_file)

echo "Checking for longRead2PairedReads.py from Battery/SciaraTools"
echo "Found: " $(which longRead2PairedReads.py)

echo "Checking for fastq2faqual.py from Battery/SciaraTools"
echo "Found: " $(which fastq2faqual.py)

echo "Checking for fasta_name_changer.py from Battery/SciaraTools"
echo "Found: " $(which fasta_name_changer.py)

