cat ../sim_pb/pacbio.s.fastq ../sim_ont/ont.s.fastq > pbont.fastq
fastq2faqual.py -f pbont.fastq --fa
