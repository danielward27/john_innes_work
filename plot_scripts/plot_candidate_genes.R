library(tidyverse)
library(viridis)

# Script follows on from identify_candidate_genes.R
collapsed_df3 = read.csv("collapsed_df_100_grouped_stats.csv")

# Plot histograms of phenotype data for accessions with and without Yr28
renseq_hits_to_accession <- function(renseq_hits){
  
  hit_accession_numbers <- renseq_hits %>%
    str_split(", ") %>%
    unlist() %>%
    strsplit("_") %>%
    sapply("[[", 2) %>%
    unique()
  
  accessions <- paste("BW", hit_accession_numbers, sep = "_")
  return(accessions)
}

Yr28_accessions <- renseq_hits_to_accession(collapsed_df3[2, "renseq_hits"])



with_yr28_phenotype <- phenotype_data[phenotype_data$accession %in% Yr28_accessions,]
without_yr28_phenotype <- phenotype_data[!(phenotype_data$accession %in% Yr28_accessions),]

# Both plots overlayed

p <- ggplot() +
  facet_wrap(~isolate) +
  scale_x_continuous(breaks = seq(1,9,1)) +
  geom_histogram(aes(mean_phenotype), with_yr28_phenotype, binwidth = 1, alpha = 0.7, fill = viridis(5)[2]) +
  geom_histogram(aes(mean_phenotype), without_yr28_phenotype, binwidth = 1, alpha = 0.7, fill = viridis(5)[4]) +
  xlab("Mean phenotype") +
  ylab("Count") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16), axis.title = element_text(size = 16), axis.text = element_text(size = 14))

p
ggsave("Phenotype_hist_by_Yr28.png", plot = p, height = 5.625, width = 10, units = "in", dpi = 300)

# Are there any resistant accessions without Yr28???
without_yr28_mean_phenotype <- without_yr28_phenotype %>%
  group_by(accession) %>%
  summarise(mean = mean(mean_phenotype))

p2 <- ggplot() +
  xlab("Mean phenotype") +
  ggtitle("Phenotype histogram excluding Yr28") +
  scale_x_continuous(breaks = seq(1,9,1), limits = c(1,9)) +
  geom_histogram(mapping = aes(mean), without_yr28_mean_phenotype, binwidth = 1) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), title = element_text(size = 14), plot.title = element_text(hjust = 0.5))

p2  
ggsave("Phenotype_hist_without_Yr28.png", plot = p2, height = 5.625, width = 10, units = "in", dpi = 300)

# Does not look like there are any resistant accessions (to all isolates) without Yr28

# Do any of the candidates explain the residual resistance (not explained by Yr28) in w057?

w057_without_yr28 <- without_yr28_phenotype %>%     # Filter to w057 data
  filter(isolate == "w057")

candidate_df <- collapsed_df3 %>%
  mutate(candidate_gene = c(1:6,9,10,7,8)) %>%
  filter(Yr28_perc_id < 99)

# Only candidate gene 4 out the 6 is found in accessions without Yr28,
# can this explain the residual resistance to w057 we observed?
accessions_with_4 <- renseq_hits_to_accession(candidate_df[candidate_df$candidate_gene == 4, "renseq_hits"])

with_4_without_Yr28 <- w057_without_yr28[w057_without_yr28$accession %in% accessions_with_4,] %>%
  mutate("candidate_4" = "With candidate 4")

without_4_without_Yr28 <- w057_without_yr28[!(w057_without_yr28$accession %in% accessions_with_4),] %>%
  mutate("candidate_4" = "Without candidate 4")

candidate_4_df <- bind_rows(with_4_without_Yr28, without_4_without_Yr28)  

ggplot(candidate_4_df) +
  xlab("Mean phenotype") +
  facet_wrap(~candidate_4) +
  ggtitle("w057 phenotype seperated by presence of candidate 4") +
  scale_x_continuous(breaks = seq(1,9,1)) +
  geom_histogram(aes(mean_phenotype), candidate_4_df, binwidth = 1.5) +
  xlab("Mean phenotype") +
  ylab("Count") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

# Only one gene (Yr28) appears to give resistance. I will rerun the analysis without accessions containing Yr28, to see if I can find a gene masked by Yr28

# Used to make a list of accessions to rerun the analysis exlcuding Yr28
not_Yr28_accessions <- phenotype_data[!(phenotype_data$accession %in% Yr28_accessions),] %>%
  select(accession) %>%
  unique() %>%
  as_tibble()
write_delim(not_Yr28_accessions, "//jic-hpc-data/group-scratch/Brande-Wulff/ward/agrenseq_without_yr28/accession_list_w057.txt",
            delim = "\n", col_names = FALSE)
