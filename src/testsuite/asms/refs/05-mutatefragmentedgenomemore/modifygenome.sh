#!/bin/bash

python2.7 -c "

from Bio import SeqIO
#f=05-ecoli.mutate-frag-mutate.fasta
f='../04-fragmentedmutatedecoligenome/04-ecoli-fragmented-mutated-genome.fasta'

fas = []
for fa in SeqIO.parse(f, 'fasta'):
  fas.append(fa)

newfas = []
while fas:
  m1 = fas.pop(0)
  m2 = fas.pop(-1)
  m1.id = m1.id + '_' +  m2.id
  m1.description = ''
  m1.seq = m1.seq + m2.seq
  newfas.append(m1)


#if anything is left, gather it
newfas += fas

out = open('temp.fasta', 'w')
for fa in newfas:
  SeqIO.write(fa, out, 'fasta')  
"




IN='temp.fasta'
OUT=05-ecoli.mutate-frag-mutate
if [ ! -f $IN.fai ]; then samtools faidx $IN; fi

../../../third_party/SVsim/SVsim -i sv.commands -r $IN -o $OUT -W


rm temp.fasta*
