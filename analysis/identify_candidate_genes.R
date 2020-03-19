library(seqinr)
library(tidyverse)

# Read in phenotype data
phenotype_files <- list.files("C:/Users/ward/Documents/AgRenSeq/for_loop_figures/phenotype_files", pattern = "phenotype_", full.names = TRUE)

read_in_phentype_file <- function(phenotype_file){
  #Get isolate to add as ID
  split_file_name <- phenotype_file %>%
    str_split("_") %>%
    unlist()
  isolate <- strsplit(split_file_name[5], ".txt")
  
  read_delim(phenotype_file, "\t", col_names = FALSE) %>%
    mutate(isolate = as.character(isolate))
}

phenotype_data <- phenotype_files %>%
  map(read_in_phentype_file) %>%
  reduce(rbind)

# Convert to tidy (long) data format
phenotype_data <- phenotype_data %>%
  gather(repeats, phenotype, 2:5) %>%
  select(accession = X1, isolate, phenotype) %>% # I won't keep which repeat is which as I'll reduce this down to means anyway
  arrange(accession)

# Taking the means
phenotype_data <- phenotype_data %>%
  group_by_at(c("accession", "isolate")) %>%
  summarize(mean_phenotype = mean(phenotype, na.rm = TRUE))

#################################################################################

# Read in AgRenSeq data

# List AgRenSeq output filenames
filenames <- list.files(pattern = "output_")

read_output_file <- function(file_name){
  #Get isolate from filename
  split_file_name <- basename(file_name %>% str_split('_',3, simplify = TRUE))
  isolate <- unlist(split_file_name)[2]
  
  # Get accession name
  start_index <- nchar(file_name) -11
  end_index <- nchar(file_name) - 4
  accession <- substring(file_name, start_index, end_index)
  
  df <- read_delim(file_name, delim = "\t", col_names = c("contig", "contig_no", "score", "kmers")) %>%
    select(-contig_no) %>%                            # Don't need the contig number (was just used for plotting)
    filter(score > 25) %>%                            # We will just look at contigs with a score >25
    mutate(accession = accession, isolate = isolate)
  return(df)
}

# Initiate df to hold all the contigs over 25
agrenseq_df <- data.frame()

# Now in a for loop add the data
for (i in 1:length(filenames)){
  df_to_add <- read_output_file(filenames[i])
  agrenseq_df <- bind_rows(agrenseq_df, df_to_add)
}

# Filter df to information important for analysis

# We will not keep all the kmer information for simplicity
# (can check graphs to visualise this), but will keep the top score and kmer number.

agrenseq_df <- agrenseq_df %>%
  mutate(ID = paste(accession, contig, sep = "_")) %>%
  group_by(ID, isolate) %>%
  filter(score == max(score)) %>%  # Gets top score for each ID and isolate
  group_by(ID) %>%
  filter(n() == 8) %>%             # We will just look at contigs that were significant across isolates (intersect)
  filter(score == max(score)) %>%  # Top score for each ID
  ungroup()
  

# Add NLR information

#Initiate collumns
agrenseq_df$complete <- ""
agrenseq_df$motif_list <- ""
agrenseq_df$nlr_length <- 0

# Read in NLR parser data and append to columns
for (i in 1:nrow(agrenseq_df)){
  accession_name <- agrenseq_df$accession[i]
  contig_name <- agrenseq_df$contig[i]
  nlr_filepath <- sprintf("//jic-hpc-data/group-scratch/Brande-Wulff/ward/genotype_data/nlrparser/nlr-parser_out_%s/nlr.txt", accession_name)
  
  nlr <- read_delim(nlr_filepath, delim = "\t", col_types = cols(MotifList = col_character())) %>%  # Have to specify MotifList as character or it just gives a long number
    filter(SequenceName == contig_name) %>%                                    # Retrieve row with contig
    mutate(nlr_length = end - start) %>%                             # Get nlr length
    select(contig = SequenceName, complete, MotifList, nlr_length)        
  
  agrenseq_df$complete[i] <- as.character(select(nlr, complete))
  agrenseq_df$motif_list[i] <- as.character(select(nlr, MotifList))
  agrenseq_df$nlr_length[i] <- as.numeric(select(nlr, nlr_length))
}

########################################################################## skip to next hash line if re-running

# Use seqinr to get contig length and extract fasta files. This is quite slow.
# If it becomes inconveniant, perhaps better to use SAMtools faidx

agrenseq_df$contig_length <- 0

for (i in 1:nrow(agrenseq_df)){
  contig <- agrenseq_df$contig[i]
  accession <- agrenseq_df$accession[i]
  ID <- paste(accession, contig, sep = "_")
  
  # Note this loop rereads the assembly each time before finding the contig, so could be worth speeding up with an if statement?
  fasta <- read.fasta(sprintf("C:/Users/ward/PycharmProjects/prepare_fastas/edited_assemblies/edited_assembly_%s.clc.fasta", accession))[paste(contig, "_", accession, sep = "", collapse = "")]
  fasta_length <- getLength(fasta)
  agrenseq_df$contig_length[agrenseq_df$ID == ID] <- fasta_length   # add contig length to df
  
  write.fasta(fasta, paste(accession, "_", contig, sep = "", collapse = ""),
              sprintf("C:/Users/ward/Documents/AgRenSeq/Analysis for presentation/significant_sequences/%s.fasta", paste(accession, "_", contig, sep = "", collapse = "")))
  print(i)
  
}


###########################################################################


# Make a significant contig list (will be used to make bash scipts)
new_sig_contigs_list <- as.data.frame(paste(agrenseq_df$accession, agrenseq_df$contig, sep = "_"))
write_delim(new_sig_contigs_list, "//jic-hpc-data/group-scratch/Brande-Wulff/ward/agrenseq_analysis/blast_sig_contigs/new_sig_contigs_list.txt",
            delim = "\t", col_names = FALSE)


# Using a python script I made job files to blast every contig in the dataframe against:
# chinese spring (CS), Yr28, BW_01111 WGS data, NWVB01 assemblies and against the other concatenated RenSeq assemblies.


# We will collate the key information from the blast files using the functions below

# RenSeq blast sseqids are formatted contig_accession, function switches this format to accession_contig
accession_then_contig <- function(contig_accession){  
  split_info <- contig_accession %>%
    as.character() %>%
    strsplit(split = "_") %>%
    unlist()
  accession_contig <- paste(split_info[c(4,5,1:3)], collapse = "_")
  
  return(accession_contig)
} 


# Function reads in the blast data (sseqid, perc_id, perc_coverage) for a particular ID and kmer_length.
# Also adds the wgs contig length. Returns a single row df.
read_wgs_blast <- function(ID){                     
  contig_length <- agrenseq_df$contig_length[agrenseq_df$ID == ID]
  
  wgs_blast <- read_delim(sprintf("//jic-hpc-data/group-scratch/Brande-Wulff/ward/agrenseq_analysis/blast_sig_contigs/blast_new_sig_wgs/blast_%s_wgs.txt", ID),
                          "\t", col_names = FALSE)[1,2:4] %>%
    rename(wgs_seqid = X2, wgs_perc_ID = X3, wgs_align_length = X4)
  
  return(wgs_blast)
}


# Reads key renseq blast info, calculates the mean phenotype of accessions that share the ID contig,
# and lists those accessions. Returns single row df.
read_renseq_blast <- function(ID, min_perc_id = 100, min_query_cov = 80){  # Defualt takes hits with 100% identity over 80% of the query length
  contig_length <- agrenseq_df$contig_length[agrenseq_df$ID == ID]
  renseq_blast <- read_delim(sprintf("//jic-hpc-data/group-scratch/Brande-Wulff/ward/agrenseq_analysis/blast_sig_contigs/blast_new_sig_contig_combined_tauschii_output/blast_new_%s_combined_tauschii.txt", ID), "\t", col_names = FALSE) %>%
    mutate(renseq_q_perc_coverage = (X4/contig_length)*100) %>%
    select(renseq_contig = X2, renseq_perc_ID = X3, renseq_q_perc_coverage) %>%
    filter(renseq_perc_ID >= min_perc_id, renseq_q_perc_coverage>min_query_cov) %>%
    rowwise() %>%
    mutate(renseq_contig = accession_then_contig(renseq_contig))

  number_with_contig <- nrow(renseq_blast)
  accessions_with_contig <- str_sub(renseq_blast$renseq_contig , 1, 8)
  phenotype_with_contig <- phenotype_data[phenotype_data$accession %in% accessions_with_contig, ]
  mean_phen_with_contig <- mean(phenotype_with_contig$mean_phenotype)
  
  IDs_with_contig_csv <- renseq_blast$renseq_contig %>%
    str_sort() %>%
    paste(collapse = ", ")
  
  renseq_blast_row <- data.frame(renseq_hits = IDs_with_contig_csv, number_with_contig = number_with_contig, mean_phen_with_contig = mean_phen_with_contig)
  return(renseq_blast_row)
}




# For each accession, read in the useful blast data using functions above, and append into data frame (one row per ID)
blast_outputs_df <- tibble() 

for (i in 1:nrow(agrenseq_df)) {
  ID <- agrenseq_df$ID[i]
  contig_length <- agrenseq_df$contig_length[agrenseq_df$ID == ID]
  
  #read and filter blast for one ID
  CS_blast <- read_delim(sprintf("//jic-hpc-data/group-scratch/Brande-Wulff/ward/agrenseq_analysis/blast_sig_contigs/blast_new_sig_contig_chinese_spring_output/blast_new_%s_chinese_spring.txt", ID), "\t", col_names = FALSE)[1,] %>% 
    select(ID = X1, chin_spring = X2, chin_spring_sstart = X9, chin_spring_send = X10) # Take useful collumns and name
  
  renseq_blast <- read_renseq_blast(ID, 99, 80)
  
  Yr28_blast <- read_delim(sprintf("//jic-hpc-data/group-scratch/Brande-Wulff/ward/agrenseq_analysis/blast_sig_contigs/blast_new_sig_contig_Yr28_output/blast_new_%s_Yr28.txt", ID), "\t", col_names = FALSE)[1,] %>%
    select(Yr28_perc_id = X3, Yr28_length = X4)
  
  NWVB01_blast <- read_delim(sprintf("//jic-hpc-data/group-scratch/Brande-Wulff/ward/agrenseq_analysis/blast_sig_contigs/blast_new_sig_contig_NWVB01_output/blast_new_%s_NWVB01.txt", ID),
                             "\t", col_names = FALSE)[1,] %>%
    select(NWVB01_seqid = X2, NWVB01_start = X9, NWVB01_end = X10)
  
  wgs_blast <- read_wgs_blast(ID)
  
  # collate the blast dfs for the ID
  blast_output_row <- tibble()
  blast_output_row <- bind_cols(CS_blast, renseq_blast, Yr28_blast, wgs_blast, NWVB01_blast) %>%
    mutate(renseq_hits = as.character(renseq_hits)) # Default was factor
  
  blast_outputs_df <- bind_rows(blast_outputs_df, blast_output_row)
}



# Bind the blast outputs with the agrenseq_df

agrenseq_df <- left_join(agrenseq_df, blast_outputs_df, by = "ID")

# Saves newline seperated wgs seqids - used to make a bash script to fetch the sequences for extending the sequence
agrenseq_df %>%
  select(wgs_seqid) %>%
  unique() %>%
  write_delim("//jic-hpc-data/group-scratch/Brande-Wulff/ward/agrenseq_analysis/BW_01111_fetch_sequences/wgs_sig_contigs.txt", delim = "\n", col_names = FALSE)

# Add phenotype data for individual accession
mean_across_accessions_phenotype_data <- phenotype_data %>%
  group_by(accession) %>%
  summarise(mean(mean_phenotype)) %>%
  rename_at(vars("mean(mean_phenotype)"), ~"mean_phenotype")

agrenseq_df <- merge(agrenseq_df, mean_across_accessions_phenotype_data, by = "accession")
  
colnames(agrenseq_df)

# Reorder dataframe
agrenseq_df <- agrenseq_df[,c(6,25,1,2:4,7:24)]

colnames(agrenseq_df)

# We now have a good df to work with
write_csv(agrenseq_df, "agrenseq_df.csv")

##########################################################################

# agrenseq_df <- read_csv("agrenseq_df.csv", col_types = cols(motif_list = col_character()))



# Function retrieves rows of agrenseq_df that match the ID contig (i.e. are in the renseq_hits collumn)
# Can be useful to find the longest contig from a resistant accession for example

ID_to_renseq_hits_df <- function(ID){
  renseq_hits <- agrenseq_df$renseq_hits[agrenseq_df$ID == ID] %>%
    str_split(pattern = ", ") %>%
    unlist()
  
  rows_index_with_hits <- vector()
  for (i in 1:length(renseq_hits)){
    rows_index_with_hits <- c(rows_index_with_hits, str_which(renseq_hits[i], agrenseq_df$ID)) 
  }
  
  hits_df <- agrenseq_df[rows_index_with_hits, ] %>%
    arrange(desc(contig_length), mean_phenotype)
  
  return(hits_df)
}


################################################################

# Reduce dataframe down to candidates:

# Find rows with overlapping hits

overlapping_hits <- function(dataframe, ID) {     # Returns dataframe of overlapping hits for an ID
  overlapping_df <- data.frame()
  renseq_hit_vector <- dataframe$renseq_hits[dataframe$ID == ID] %>%
    str_split(pattern = ", ") %>%
    unlist()
  
  for (i in 1:nrow(dataframe)){
    
    check_match_vector <- dataframe$renseq_hits[i] %>%
      str_split(pattern = ", ") %>%
      unlist()
    
    intersect_length <- length(intersect(renseq_hit_vector, check_match_vector))
    if (intersect_length > 0){
      overlapping_df <- bind_rows(overlapping_df, dataframe[i,])
    }
  }
  return(overlapping_df)
}   


# Takes rows with overlapping renseq hits and collapses them
# (taking most resistant, but with an inclusive renseq_hits collumn, mean phenotype and number with contig)
collapse_df <- function(dataframe){  
  collapsed_df <- data.frame()
  
  # For all IDs
  for (i in 1:nrow(dataframe)){
    overlapping_df <- overlapping_hits(dataframe, dataframe$ID[i])
    
    # Reduce down to one row
    
    shared_renseq_hits <- overlapping_df$renseq_hits %>%
      str_split(", ") %>%
      unlist() %>%
      unique() %>%
      sort()
    
    split_shared_renseq_hits <- shared_renseq_hits %>%
      str_split("_") %>%
      sapply("[[",2)
    
    shared_accessions <- paste("BW", split_shared_renseq_hits, sep = "_")
    phenotype_with_contig <- phenotype_data[phenotype_data$accession %in% shared_accessions, ]
    shared_mean_phen_with_contig <- mean(phenotype_with_contig$mean_phenotype)
  
    shared_number_with_contig <- length(unique(shared_accessions))
    
    shared_renseq_hits <- paste(shared_renseq_hits, collapse = ", ")
    
    most_resistant <- overlapping_df %>%
      filter(mean_phenotype == min(mean_phenotype)) %>%
      filter(contig_length == max(contig_length))
    
    single_row <- most_resistant[1,] %>%
      mutate(renseq_hits = shared_renseq_hits,
             number_with_contig = shared_number_with_contig,
             mean_phen_with_contig = shared_mean_phen_with_contig)
    
    
    collapsed_df <- bind_rows(collapsed_df, single_row)  
  }
  
  collapsed_df <- distinct(collapsed_df)
  return(collapsed_df)
}


collapsed_df <- collapse_df(agrenseq_df)
collapsed_df2 <- collapse_df(collapsed_df) 
collapsed_df3 <- collapse_df(collapsed_df2) # Reapplying means that if A matches B, which matches C, but A doesn't reach the threshold to match C, they will still be grouped as a single group.

write_csv(collapsed_df3, "collapsed_df_100_grouped_stats.csv") # These are the candidate genes

#################################################################

# Summary statistics:

agrenseq_df %>%
  select(complete) %>%
  group_by(complete) %>%
  tally()

agrenseq_df %>%
  select(contig_length) %>%
  summarise(mean = mean(contig_length), sd = sd(contig_length))


##################################################################

