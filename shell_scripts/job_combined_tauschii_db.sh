#!/bin/bash
#SBATCH --job-name=combined_tauschii_db
#SBATCH -o combined_tauschii_db.out
#SBATCH -e combined_tauschii_db.err
#SBATCH --mem 18gb
#SBATCH -p nbi-short   # Partition (queue equivalent)

export SBATCH_PARTITION=nbi-short
source blast+-2.7.1

makeblastdb -in /jic/scratch/groups/Brande-Wulff/ward/genotype_data/assemblies/combined_tauschii.clc.fasta -dbtype nucl -out combined_tauschii
