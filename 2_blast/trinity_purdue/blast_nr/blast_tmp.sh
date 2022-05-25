#!/bin/bash
#SBATCH -N 1
#SBATCH -n 128
#SBATCH --exclusive
#SBATCH -A standby
#SBATCH --time=4:00:00
#SBATCH --job-name=blastx
#SBATCH --error=%x-%J-%u.err
#SBATCH --output=%x-%J-%u.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lee3617@purdue.edu

# Prerequisite modules
ml utilities hyper-shell 

# Make a copy of database and input in /tmp.
# TODO: possibly add a check for sufficient /tmp free space later.
# Or just error out if the copy fails.

FASTA_FILE=$1
echo $FASTA_FILE >> "test_fasta.txt"
myseqID="$(basename -- $FASTA_FILE .fa)"
tempdir=/tmp/$USER/blastx;mkdir -p $tempdir/{input,database} || exit 1
cd $tempdir

cp -R $RCAC_SCRATCH/KW/blast_nr/database/. $tempdir/database/ || exit 1

cp $RCAC_SCRATCH/KW/blast_nr/$FASTA_FILE $tempdir/ || exit 1

#Split the fasta file
awk -F "|" '/^>/ {close(F); ID=$1; gsub("^>", "", ID); F="input/"ID".fasta"} {print >> F}' $FASTA_FILE || exit 1
for FILE in input/*.fasta
 do
    ID="$(basename -- $FILE .fasta)"
    Output=$ID"_blastxout"
    echo "/group/bioinfo/apps/apps/blast-2.12.0+/bin/blastx  -num_threads 1 -evalue 1e-3 -query "/tmp/$USER/blastx/input/"$ID".fasta" -db "/tmp/$USER/blastx/database/"Lophotrochozoa_model -out $RCAC_SCRATCH"/KW/blast_nr/${myseqID}_blastoutput/"$Output -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname scomname stitle\" && echo $ID success" >> ${tempdir}/blastTaskfile.txt
 done || exit 1
~                 

# And go!
mkdir -p $RCAC_SCRATCH/KW/blast_nr/${myseqID}_blastoutput || exit 1

hyper-shell cluster ${tempdir}/blastTaskfile.txt \
	-o $RCAC_SCRATCH"/KW/blast_nr/${myseqID}_Taskfile.output" \
        -f $RCAC_SCRATCH"/KW/blast_nr/${myseqID}_Taskfile.failed" \
	-N 128 # -N 128 here represent 128 jobs/commands run at the same time 
