#!/bin/bash

python2.7 -c "
from Bio import SeqIO
import numpy as np
import random

f = '../01-ecoligenome/01-ecoli_k12_mg1655_NC_000913.3.fasta'
fa = SeqIO.read(f,'fasta')

snprate=1/100.0
orig = fa.id
fa.description = ''
newseq=''
for i in range(0,len(fa)):
  b = fa.seq[i]
  if np.random.binomial(1,snprate):
    newb = b
    while newb == b:
        newb = 'ACGT'[random.randint(0,3)]
    newseq += newb
  else:
    newseq += b 
fa.seq = newseq

with open('02-ecoli-snpd-genome.fasta','w') as out:
  out.write('>ecoli_snp\n'+newseq+'\n')
"
