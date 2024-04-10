# Load ggplot2
library(ggplot2)

# Assuming your data is in a file named 'data.txt'
data <- read.table('/home/guarracino/Desktop/pggb-paper/hsapiens90.chr6.masked_p98.s20000.n90.k47_1.degree-depth.chm13.10kbp.tsv', header = T, sep = '\t', comment.char = "?")

options(scipen=10000)

# Plot
p <- ggplot(data, aes(x = start, y = value, color = smoothxg, group = interaction(smoothxg, metric))) +
  geom_line() +
  facet_grid(metric ~ ., scales = "free") + 
  #theme_minimal() +
  theme(
    text = element_text(size = 24), # Default text size
    axis.title = element_text(size = 16), # Axis title text size
    axis.text = element_text(size = 16), # Axis text size
    plot.title = element_text(size = 20, hjust = 0.5), # Center the title
    legend.title = element_text(size = 16), # Legend title text size
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 16)
  ) +
  scale_y_log10() +
  labs(title = "Human chromosome 6 pangenome graph",
       x = "Position",
       y = "Value",
       color = "SMOOTHXG")
#+ xlim(29595119, 32911317)
#+ xlim(58227903, 62704682)
#+ 
path_image <- "/home/guarracino/Desktop/plot_degree_depth_log_scale.pdf"
ggsave(path_image, p, width = 13, height = 7)

