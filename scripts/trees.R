#!/usr/bin/env Rscript

library(tidyverse)
library(dendextend) # for tanglegram
library(ape) # for phylogenetic analysis
library(ggtree) # for phylogenetic tree plotting and annotation

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if we have exactly 3 arguments
if(length(args) != 3) {
  stop("Usage: Rscript make_tree.R path_dist_tsv plot_title output_file_path", call. = FALSE)
}

# Assign arguments
path_dist_tsv <- args[1]
plot_title <- args[2]
output_file_path <- args[3]

# Read matrices
jaccard_dist_df <- read_tsv(path_dist_tsv) %>%
  arrange(group.a, group.b) %>%
  select(group.a, group.b, jaccard.distance) %>%
  pivot_wider(names_from = group.b, values_from = jaccard.distance) %>%
  column_to_rownames(var = "group.a")

#euclidean_dist_df <- read_tsv(path_dist_tsv) %>%
#  arrange(group.a, group.b) %>%
#  select(group.a, group.b, euclidean.distance) %>%
#  pivot_wider(names_from = group.b, values_from = euclidean.distance) %>%
#  column_to_rownames(var = "group.a")

jaccard_hc <- as.dist(jaccard_dist_df) %>% hclust()
#euclidean_hc <- as.dist(euclidean_dist_df) %>% hclust()

# Plot
pdf(output_file_path) # Start PDF device, replace with svg(output_file_path) for SVG format
plot(
  jaccard_hc,
  # Label at same height
  hang = -1,
  main = plot_title,
  xlab = 'Haplotype',
  ylab = 'Jaccard distance',
  sub = '',
  cex = 1.8,       # Adjusts the size of points/text in the plot
  cex.lab = 1.6,   # Adjusts the size of x and y labels
  cex.axis = 1.6,  # Adjusts the size of axis text
  cex.main = 1.6,  # Adjusts the size of the main title
  cex.sub = 1.6,    # Adjusts the size of the subtitle
  lwd = 2  # Increases the width of the branch lines
)
dev.off() # Close the device


if (FALSE) {
  #library(RColorBrewer)
  #library(ggrepel)
  
  tree <- nj(as.dist(jaccard_dist_df)) # compute the tree using the nj function from the ape package
  ggtree(tree) # plot the tree using the ggtree function from the ggtree package
  
  ggtree(
    tree,
    ladderize=T,
    
    # branch.length
    branch.length="branch.length",
    
    #layout="daylight"
  ) + 
    #geom_point(aes(shape=isTip, color=isTip), size=3) + 
    geom_nodepoint(color="#b5e521", alpha=1/4, size=6) + 
    #geom_tippoint(color="#FDAC4F", shape=8, size=3) +
    geom_tiplab(size=3, color="black") +
    geom_rootedge(rootedge = 0.5) + ggplot2::xlim(0, 0.8)
  
  
  jaccard_dend <- as.dist(jaccard_dist_df) %>% hclust() %>% as.dendrogram()
  euclidean_dend <- as.dist(euclidean_dist_df) %>% hclust() %>% as.dendrogram()
  
  jaccard_dend %>% plot
  euclidean_dend %>% plot
  
  dl <- dendlist(jaccard_dend, euclidean_dend)
  tanglegram(dl, sort = TRUE, common_subtrees_color_lines = FALSE, highlight_distinct_edges  = FALSE, highlight_branches_lwd = FALSE)
  
  
  
  
  # Convert hclust into a dendrogram and plot
  
  # Define nodePar
  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                  cex = 0.7, col = "black")
  
  # Default plot
  plot(
    as.dendrogram(jaccard_hc),
    type = "rectangle",
    ylab = "Jaccard distance",
    horiz = F,
    #nodePar = nodePar,
    margin=0.2
  )
  
  
  # install.packages("ape")
  library("ape")
  # Default plot
  plot(
    as.phylo(jaccard_hc),
    cex = 0.9,
    label.offset = 0.005
  )
  
  library("ggplot2")
  library("ggdendro") # install.packages("ggdendro")
  ggdendrogram(
    jaccard_hc, rotate = F, theme_dendro = T)
  
  
  
  
  
  
  
  
  
  
  
  
  
  library(tidyverse)
  library(dendextend) # for tanglegram
  library(ape) # for phylogenetic analysis
  library(ggtree) # for phylogenetic tree plotting and annotation
  
  #library(RColorBrewer)
  #library(ggrepel)
  
  path_dist1_tsv <- '~/Desktop/scerevisiae142.fa.gz.2cc7813.e34d4cd.d05ee65.smooth.final.dist.tsv'
  path_dist2_tsv <- '~/Desktop/scerevisiae142.fa.gz.f4dcc26.e34d4cd.d05ee65.smooth.final.dist.tsv'
  
  # Read matrices
  jaccard_dist1_df <- read_tsv(path_dist1_tsv) %>%
    arrange(group.a, group.b) %>%
    mutate(jaccard=1-jaccard) %>%
    select(group.a, group.b, jaccard) %>%
    pivot_wider(names_from = group.b, values_from = jaccard) %>%
    column_to_rownames(var = "group.a")
  euclidean_dist1_df <- read_tsv(path_dist1_tsv) %>%
    arrange(group.a, group.b) %>%
    select(group.a, group.b, euclidean) %>%
    pivot_wider(names_from = group.b, values_from = euclidean) %>%
    column_to_rownames(var = "group.a")
  
  jaccard_dist2_df <- read_tsv(path_dist2_tsv) %>%
    arrange(group.a, group.b) %>%
    mutate(jaccard=1-jaccard) %>%
    select(group.a, group.b, jaccard) %>%
    pivot_wider(names_from = group.b, values_from = jaccard) %>%
    column_to_rownames(var = "group.a")
  euclidean_dist2_df <- read_tsv(path_dist2_tsv) %>%
    arrange(group.a, group.b) %>%
    select(group.a, group.b, euclidean) %>%
    pivot_wider(names_from = group.b, values_from = euclidean) %>%
    column_to_rownames(var = "group.a")
  
  
  
  tree <- nj(as.dist(jaccard_dist_df)) # compute the tree using the nj function from the ape package
  ggtree(tree) # plot the tree using the ggtree function from the ggtree package
  
  ggtree(
    tree,
    ladderize=T,
    
    # branch.length
    branch.length="none",
    
    #layout="daylight"
  ) + 
    #geom_point(aes(shape=isTip, color=isTip), size=3) + 
    geom_nodepoint(color="#b5e521", alpha=1/4, size=6) + 
    #geom_tippoint(color="#FDAC4F", shape=8, size=3) +
    geom_tiplab(size=3, color="black") +
    geom_rootedge(rootedge = 0.5) + ggplot2::xlim(-0.5, 30)
  
  tree2 <- as.dist(jaccard_dist_df) %>% hclust()
  ggtree(
    tree2,
    ladderize=T,
    
    # branch.length
    branch.length="branch.length",
    
    #layout="daylight"
  ) + 
    #geom_point(aes(shape=isTip, color=isTip), size=3) + 
    geom_nodepoint(color="#b5e521", alpha=1/4, size=6) + 
    #geom_tippoint(color="#FDAC4F", shape=8, size=3) +
    geom_tiplab(size=3, color="black") +
    geom_rootedge(rootedge = 0.5) + ggplot2::xlim(0, 0.17)
  
  
  
  jaccard_dend <- as.dist(jaccard_dist_df) %>% hclust() %>% as.dendrogram()
  euclidean_dend <- as.dist(euclidean_dist_df) %>% hclust() %>% as.dendrogram()
  
  jaccard_dend %>% plot 
  euclidean_dend %>% plot
  
  dl <- dendlist(jaccard_dend, euclidean_dend)
  tanglegram(dl, sort = TRUE, common_subtrees_color_lines = FALSE, highlight_distinct_edges  = FALSE, highlight_branches_lwd = FALSE)
  
  
  jaccard_dend1 <- as.dist(jaccard_dist1_df) %>% hclust() %>% as.dendrogram()
  jaccard_dend2 <- as.dist(jaccard_dist2_df) %>% hclust() %>% as.dendrogram()
  dl12 <- dendlist(jaccard_dend1, jaccard_dend2)
  tanglegram(dl12, sort = TRUE, common_subtrees_color_lines = FALSE, highlight_distinct_edges  = FALSE, highlight_branches_lwd = FALSE)
}
