#!/bin/bash
#SBATCH --job-name=makeblastdb_BW_01111
#SBATCH -p nbi-short
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o job_output_makeblastdb_BW_01111.out
#SBATCH -e job_error_makeblastdb_BW_01111.err
#SBATCH --mem 8gb

export SBATCH_PARTITION=nbi-short

source blast+-2.7.1
makeblastdb -in BW_01111.unitigs.fa -dbtype nucl -out BW_01111.unitigs.fa
