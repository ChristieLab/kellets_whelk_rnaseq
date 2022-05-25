#!/bin/bash

# split the fasta into three smaller fasta files containing 26000 in the first two fastas
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%26000==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < notfoundbyswiss.fasta || exit 1

# submit slurm job for each fasta file
for FILE in myseq*.fa
do
	echo $FILE 
	sbatch blast_tmp.sh $FILE
done || exit 1


