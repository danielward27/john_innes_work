#!/bin/bash
#SBATCH --job-name=index_combined_tauschii
#SBATCH -o index_combined_tauschii.out
#SBATCH -e index_combined_tauschii.err
#SBATCH --mem 10gb
#SBATCH -p nbi-short   # Partition (queue equivalent)

export SBATCH_PARTITION=nbi-short
source samtools-1.9

samtools faidx /jic/scratch/groups/Brande-Wulff/ward/genotype_data/assemblies/combined_tauschii.clc.fasta