import sys, os, re
import random
from itertools import islice

def strMUT(text, dic):
    """ Replaces keys of dic with values of dic in text. 2005-02 by Emanuel Rumpf """
    pat = "(%s)" % "|".join(map(re.escape, dic.keys()))
    return re.sub(pat, lambda m: dic[m.group()], text)


header='''#!/bin/bash
#SBATCH --job-name=blast_ACCESSION_CONTIG_Yr28
#SBATCH -p nbi-short
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o job_output_blast_ACCESSION_CONTIG_Yr28.out
#SBATCH -e job_error_blast_ACCESSION_CONTIG_Yr28.err
#SBATCH --mem 8gb

export SBATCH_PARTITION=nbi-short

source blast+-2.7.1
blastn -db /jic/scratch/groups/Brande-Wulff/ward/agrenseq_analysis/blast_sig_contigs/Yr28_database/Yr28_10kb.fasta -query /jic/scratch/groups/Brande-Wulff/ward/agrenseq_without_yr28/sig_contigs/ACCESSION_CONTIG.fasta -out /jic/scratch/groups/Brande-Wulff/ward/agrenseq_without_yr28/blast/blast_Yr28/blast_ACCESSION_CONTIG_Yr28.txt -perc_identity 1 -outfmt 6 -task blastn

'''

PARTS = []
with open("sig_contigs_list.txt", "r") as f:   
	for line in f:
		PARTS.append(line.rstrip())          # adds ACCESSION_CONTIGs to ACCESSION_CONTIG list

with open("submitter_blast_sig_Yr28.sh" ,'w') as subOutF: 
 subOutF.write('#!/bin/bash\n')
 for pita, PART in enumerate(PARTS,1):##
   outF='DAN_blast_sig' + str(pita) + '.sh'
   subOutF.write('sbatch '+ outF + '\n')
   with open(outF,'w') as out:
     newtemplate=strMUT(header, {'ACCESSION_CONTIG':  str(PART)})
     out.write(newtemplate)
  
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 