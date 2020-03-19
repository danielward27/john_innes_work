import os
from Bio import SeqIO

# get file and accession list
assembly_filenames = os.listdir("assemblies")

accession_list = []
for i in range(len(assembly_filenames)):
    filename_replace = assembly_filenames[i].replace('.', '_')
    filename_split = filename_replace.split('_')
    accession_name = filename_split[1] + "_" + filename_split[2]
    accession_list.append(accession_name)

# Add accession onto record ID of fasta files
for i in range(len(assembly_filenames)):
    original_assembly = "assemblies/" + assembly_filenames[i]
    edited_assembly = "edited_assemblies/edited_" + assembly_filenames[i]

    with open(original_assembly) as original, open(edited_assembly, 'w') as edited:
        records = SeqIO.parse(original_assembly, 'fasta')
        for record in records:
            record.id = record.id + "_" + accession_list[i]
            record.description = ""
            print(record.id)                            # not required but allows you to track progress
            SeqIO.write(record, edited, 'fasta')




