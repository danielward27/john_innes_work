#!/bin/bash
#SBATCH --job-name=samtools_BW_01111_index
#SBATCH -o samtools_BW_01111_index.out
#SBATCH -e samtools_BW_01111_index.err
#SBATCH --mem 4gb
#SBATCH -p nbi-short   # Partition (queue equivalent)

export SBATCH_PARTITION=nbi-short
source samtools-1.9

samtools faidx /jic/scratch/groups/Brande-Wulff/ward/genotype_data/assemblies/assembly_BW_01111.clc.fasta