#Script will plot NLR distribution across the chromosomes of Aegilops sharonensis

# Install packages if not installed
library(tidyverse)
library(viridis)

# Read in chromosome length
df = read_csv("Table1_AS_1644_Assembly_chromosome_length.csv", col_names = TRUE)

# Change format of length from Mb to bp (e.g. 700 Mb to 700000000)
reformat_length <- function(length){
  split_length <- length %>%
    str_split(" ") %>%
    unlist()
  
  lengthbp <- as.numeric(split_length[1])*10^6 # convert to bp
  return(lengthbp)
}

df <- df %>%
  rowwise() %>%
  mutate(length_bp = reformat_length(length)) %>%      # Reformat length
  select(chr, length_bp)

# Width parameters
chr_width <- 7*10^7
x1_positions <- seq(0, 18*chr_width, 3*chr_width)

# x coordinates of the chromosomes
df$x1 <- x1_positions                             
df$x2 <- x1_positions + chr_width
df$y2 <- max(df$length_bp)
df$y1 <- df$y2 - df$length_bp

# Now read in NLR location
df_pos <- read_csv("1S_7S_complete_nlr_coordinates.csv") %>%
  mutate(chr = str_sub(chr, 4)) # remove "chr" preceding chromosome ID

df_pos <- left_join(df_pos, df, by = "chr") %>%    # Join the dataframes
  mutate(position = max(length_bp)-position)       # Make the positions top down

# Uses viridis colour package - so should print fine in black and white and is colourblind friendly
fill <- viridis(9, option = "D")[3]
line_colour <- viridis(2, option = "C")[2]

p <- ggplot(df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, label = chr)) +
  geom_rect(fill = fill) + # Can add colour = "black" to add border
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank()) +  # Remove axis and gridlines etc.
  geom_text(aes(x=x1,y=y2), nudge_x = 3.4*10^7, nudge_y = 3.7*10^7, size = 5) +     # Add and adjust position and size of label
  geom_segment(data = df_pos, aes(x = x1, xend = x2, y = position, yend = position),
               colour = line_colour) + # Adds lines for NLR distribution
  geom_segment(data = df, aes(x = (max(x2)+0.8*chr_width), xend = (max(x2)+0.8*chr_width),
                              y = y2 - 10^8, yend = y2), colour = "black", size = 0.7) +
  geom_text(aes(x=max(x1), y=y2, label = "100 Mb"), nudge_x = 2.5*chr_width, nudge_y = -4.5*10^7, angle = 270, size = 4.7)
p

ggsave("rplot.png", p, dpi = 300, width = 15, height = 15, units = "cm")

# With border
p + geom_rect(colour = "black", fill = "transparent", size = 0.5)


# Zebra version
fill = "white"
line_colour = "black"

ggplot(df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, label = chr)) +
  geom_rect(fill = fill) + # Can add colour = "black" to add border
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank()) +  # Remove axis and gridlines etc.
  geom_text(aes(x=x1,y=y2), nudge_x = 3.4*10^7, nudge_y = 4*10^7) +     # Add and adjust position of label
  geom_segment(data = df_pos, aes(x = x1, xend = x2, y = position, yend = position),
               colour = line_colour) + # Adds lines for NLR distribution
  geom_segment(data = df, aes(x = (max(x2)+0.8*chr_width), xend = (max(x2)+0.8*chr_width),
                              y = y2 - 10^8, yend = y2), colour = "black", size = 0.7) +
  geom_text(aes(x=max(x1), y=y2, label = "100 Mb"), nudge_x = 2.4*chr_width, nudge_y = -4*10^7, angle = 270) +
  geom_rect(colour = line_colour, fill = "transparent")

