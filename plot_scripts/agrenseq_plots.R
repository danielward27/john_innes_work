setwd("C:/Users/ward/Documents/AgRenSeq/ggplot2_practice")
library(ggplot2)
library(dplyr)


# Read in output add collumns names
read_output <- function(output_file){
  output_df <- read.csv(output_file, sep = "\t", header = FALSE,
                        col.names = c("contigs", "x_value", "score", "kmers"))
  output_df <- mutate(output_df, kmerXscore = kmers*score)       # Adds kmer X score as indicator of significance
  output_df <- mutate(output_df, point_size = log(kmers, 2)/2)   # Adds reasonable point size for plotting
}


# Adds binary column for contigs over threshold score and kmerXscore
add_sig_contigs <- function(output_df, min_score = 30, min_kmerXscore = 0){
  output_df$significant_contigs <- "0"
  sig_score_contigs <- as.vector(unique(output_df$contigs[output_df$score>min_score])) # Get vector for contigs which exceed score threshold
  kmerXscore_sum_list <- list()                                                        # Get kmerXscore sum for each contig above score threshold
  for (contig in sig_score_contigs){
    sig_score_df <- filter(output_df, contigs == contig)
    kmerXscore_sum <- sum(sig_score_df$kmerXscore)
    kmerXscore_sum_list <- c(kmerXscore_sum_list, kmerXscore_sum)
  }
  kmerXscore_sum_list
  significant_contigs <- sig_score_contigs[which(kmerXscore_sum_list>min_kmerXscore)] # sig kmerXscore and kmer number
  for (i in 1:nrow(output_df)) {
    for (contig in significant_contigs){
      if (output_df$contigs[i] == contig){
        output_df$significant_contigs[i] = "1"
      }
    }
  }
  return(output_df)
}

# plot one graph
output_df <- read_output("without_yr28_output_w057_BW_01179.txt")  
unique(output_df$contigs[output_df$score>30])

output_df <- add_sig_contigs(output_df, 25, 0)      # to colour over threshold 25

my_contigs <- c("contig_399_1", "contig_90_1")

# to instead colour specific contigs
output_df$significant_contigs <- "0"
output_df$significant_contigs <- output_df$contigs %in% my_contigs
output_df$labels <- ""   # Add column for labels
output_df$labels[output_df$significant_contigs == TRUE] <- "5D"

ggplot(output_df, aes(x = x_value, y = score, colour = factor(significant_contigs))) +
  theme(legend.position = "none") +
  geom_point(size = output_df$point_size, alpha = 0.7) +
  labs(x = "contigs", title = "Accession = BW_01179  -  Pathogen isolate = w057") +
  geom_text(aes(x = x_value, y = max(score), label = output_df$labels), nudge_x = 50)


#  Make a background for a slide
output_df <- output_df %>%
mutate(point_size = point_size*3)

ggplot(output_df, aes(x = x_value, y = score, colour = factor(significant_contigs))) +
theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank(), panel.background = element_rect(fill = hsv(0.56,0.6,0.33))) +
  scale_color_manual(values = c("deepskyblue3","#E69F00" )) +
  geom_point(size = output_df$point_size, alpha = 0.7)
  

################ for loop plots ####################

filenames <- list.files(pattern = "output_")
filenames

#get accession list and strain list

splitfilenames <- strsplit(filenames, '.', fixed = TRUE)
splitfilenames <-unlist(splitfilenames)
names_index <- seq(1, length(splitfilenames), 2)
splitfilenames <- splitfilenames[names_index]
splitfilenames <- strsplit(splitfilenames, '_')
splitfilenames <- unlist(splitfilenames)
head(splitfilenames)

strain_vector <- splitfilenames[seq(2, length(splitfilenames), 4)]
strain_vector

accession_vector <- paste("BW_", splitfilenames[seq(4, length(splitfilenames), 4)], sep = "")
accession_vector

for (i in 1:length(filenames)){
  png(paste(c(accession_vector[i], "_", strain_vector[i], ".png"), collapse = ""), 10, 5, units = "in", res = 600)
  
  output_df <- read_output(filenames[i])  
  output_df <- add_sig_contigs(output_df, 25, 0)
  
  graph <- ggplot(output_df, aes(x = x_value, y = score, colour = factor(significant_contigs))) +
    geom_point(size = output_df$point_size, alpha = 0.7) +
    labs(x = paste(c("NLR Contigs - Accession ", accession_vector[i]), collapse = ""),
         title = paste(c("Pathogen isolate", strain_vector[i], "- Accession", accession_vector[i]), collapse = " ")) +
    theme(legend.position = "none",
          axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 11, face = "bold"),
          title = element_text(11, face = "bold"))
  print(graph)
  dev.off()
}
  


# Make faceted example graphs for presentation

PSTv04_df <- read_output("output_PSTv04_BW_01068.txt") %>%
  mutate(isolate = "PSTv04")
PSTv18_df <- read_output("output_PSTv18_BW_01068.txt") %>%
  mutate(isolate = "PSTv18")
PSTv51_df <- read_output("output_PSTv51_BW_01068.txt") %>%
  mutate(isolate = "PSTv51")
w010_df <- read_output("output_w010_BW_01068.txt") %>%
  mutate(isolate = "w010")
w011_df <- read_output("output_w011_BW_01068.txt") %>%
  mutate(isolate = "w011")
w034_df <- read_output("output_w034_BW_01068.txt") %>%
  mutate(isolate = "w034")
w056_df <- read_output("output_w056_BW_01068.txt") %>%
  mutate(isolate = "w056")
w057_df <- read_output("output_w057_BW_01068.txt") %>%
  mutate(isolate = "w057")

facet_df <- bind_rows(PSTv04_df, PSTv18_df, PSTv51_df, w010_df, w011_df, w034_df, w056_df, w057_df, ) %>%
  mutate(isolate_contig = paste(isolate, contigs, sep = "_")) # Means we can colour individually for each graph easily

facet_df$significant_contigs <- "0"

# Get vector for contigs which exceed score threshold
sig_score_contigs <- facet_df %>%
  filter(score >=25) %>%
  select(isolate_contig) %>%
  unique() %>%
  unlist() %>%
  as.vector()

for (i in 1:nrow(facet_df)) {
  for (contig in sig_score_contigs){
    if (facet_df$isolate_contig[i] == contig){
      facet_df$significant_contigs[i] = "1"
    }
  }
}


p <- ggplot(facet_df, aes(x = x_value, y = score, colour = factor(significant_contigs))) +
  geom_point(size = facet_df$point_size, alpha = 0.7) +
  theme_bw() +
  labs(title = "BW_01068", x = "BW_01068 contigs", y = "Score") +
  facet_wrap("isolate") +
  theme(legend.position = "none", title = element_text(size = 16), plot.title = element_text(hjust = 0.5), strip.text = element_text(size = 17), axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17), axis.text = element_text(size=15))

ggsave("facet_BW_01068.png", p, width = 50, height = 28.125, units = "cm", dpi = 200)


