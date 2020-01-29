#!/bin/bash

python -c "
from Bio import SeqIO
import random

f = '../../refs/ecoligenome/01-ecoli_k12_mg1655_NC_000913.3.fasta' 
fa = SeqIO.read(f,'fasta')

lenfa = len(fa)
A=500
B=1500

orig = fa.id
fa.description = ''
with open('inputasm.fasta','w') as out:
  i=1
  lastpos = 0
  pos = lastpos + random.randint(A,B)
  fa.id = orig + '_' + str(i)
  SeqIO.write(fa[lastpos:pos], out, 'fasta')
  while pos < lenfa:
    i+=1
    lastpos = pos
    pos = lastpos + random.randint(A,B)
    fa.id = orig + '_' + str(i)
    SeqIO.write(fa[lastpos:pos], out, 'fasta')
  
"
