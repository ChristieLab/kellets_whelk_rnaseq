#!/bin/bash

# split the fasta into smaller fasta files containing 20000 contigs (rcac suggested 26000 for trinity) 
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%20000==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < rnaspades_assembly_calpoly.fasta || exit 1

# submit slurm job for each fasta file
for FILE in myseq*.fa
do
	echo $FILE 
	sbatch blast_tmp.sh $FILE
done || exit 1


