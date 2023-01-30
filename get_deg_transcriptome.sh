#!/bin/bash
#SBATCH --job-name=fgrep
#SBATCH -A beagle
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lee3617@purdue.edu
#SBATCH --output slurm-%x-%J.out

cd $SLURM_SUBMIT_DIR
module purge

tail -2770 degs_padj05.MON.NAP.csv | cut -d ',' -f8 | sed -e 's/^"//' -e 's/"$//' | sed 's/^/>/' > rnaspades_degs_contignames.txt
fgrep -f rnaspades_degs_contignames.txt -A1 --no-group-separator rnaspades_assembly_calpoly_singleline.fasta > rnaspades_degs_transcriptome.fasta
