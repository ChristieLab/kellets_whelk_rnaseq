#!/bin/bash
#SBATCH --job-name=ranger_loo
#SBATCH --array=1-60
#SBATCH -A beagle
#SBATCH --mem=9G
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 4-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lee3617@purdue.edu
#SBATCH --output slurm-%x-%J.out
module purge
module load r/4.3.1

data=./all_nocross.rds
out=./ranger_loo

Rscript ranger_loo.R $data $out $SLURM_ARRAY_TASK_ID
