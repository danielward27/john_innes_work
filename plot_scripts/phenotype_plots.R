setwd("C:/Users/ward/Documents/AgRenSeq/phenotype_correlation")
library(tidyverse)
library(corrr)
library(corrplot)
library(viridis)
library(RColorBrewer)

get_means <- function(phenotype_file, isolate) {
  df <- read.table(phenotype_file, sep="\t")
  df_mean <- data.frame("accession" = df$V1,
                            mean = rowMeans(df[,2:5], na.rm = TRUE)) %>%
    mutate(isolate = isolate)
}

PSTv04 <- get_means("phenotype_PSTv04.txt", "PSTv04")
PSTv18 <- get_means("phenotype_PSTv18.txt", "PSTv18")
PSTv51 <- get_means("phenotype_PSTv51.txt", "PSTv51")
w010 <- get_means("phenotype_w010.txt", "w010")
w011 <- get_means("phenotype_w011.txt", "w011")
w034 <- get_means("phenotype_w034.txt", "w034")
w056 <- get_means("phenotype_w056.txt", "w056")
w057 <- get_means("phenotype_w057.txt", "w057")

# Bind all the dfs

df <- bind_rows(PSTv04, PSTv18, PSTv51, w010, w011, w034, w056, w057)
# Plot histograms for the phenotype data
df %>%
  group_by(accession) %>%
  summarise(mean = mean(mean)) %>%
  filter(mean < 3) %>%
  print(n = 20)

p <- ggplot(df, aes(x=mean)) +
  theme_bw() +
  geom_histogram(colour = "black", binwidth = 1) +
  facet_wrap(~ isolate) +
  xlab("Mean phenotype") +
  ylab("Accessions count") +
  scale_x_continuous(breaks = seq(0,9,1)) +
  theme(strip.text.x = element_text(size = 14), axis.text = element_text(size = 12), axis.title = (element_text(size=16)))
p

ggsave("phenotype_histogram.png", p, dpi = 300, width = 10, height = 5.625, units = "in")

correlations <- cor(means_df[2:9], use = "all.obs", method = "pearson")

#### Correlations plot

# Calulate correlations
cor_df <- spread(df, isolate, mean) %>%
  select(-accession) %>%
  cor(use='pairwise')

col <- brewer.pal(n=11, name="RdBu")
col <- c(rep(col[1], 4), col) # Shifts scale to the right as my correlations do not vary from -1 to 1

# Correlations matrix plot
corrplot(cor_df, method = "shade", addCoef.col = "black", cl.lim=c(0.6,1), tl.col = "black",
         tl.srt=50, type = "lower", is.corr = FALSE,  col = col)


ggplot(df, aes(x=mean)) +
  theme_bw() +
  scale_x_continuous(breaks = seq(1,9,1)) +
  geom_density(size = 0.8) +
  ylab("Accession density") +
  xlab("Mean phenotype score") +
  theme(axis.title = element_text(size = 23), axis.text = element_text(size = 19))

