library(tidyverse)
library(dendextend) # for tanglegram
library(ape) # for phylogenetic analysis
library(ggtree) # for phylogenetic tree plotting and annotation BiocManager::install("ggtree")

path_dist_tsv <- '~/Desktop/primates16.hsa6_p90.s10000.n16.k47.G700-900-1100.Pasm5.O0001_1.dist.tsv'

# Read matrices
jaccard_dist_df <- read_tsv(path_dist_tsv) %>%
  #filter(! group.a %in% c('mSymSyn1#1', 'mSymSyn1#2')) %>%
  #filter(! group.b %in% c('mSymSyn1#1', 'mSymSyn1#2')) %>%
  arrange(group.a, group.b) %>%
  #mutate(jaccard=1-jaccard) %>%
  select(group.a, group.b, jaccard.distance) %>%
  pivot_wider(names_from = group.b, values_from = jaccard.distance) %>%
  column_to_rownames(var = "group.a")
euclidean_dist_df <- read_tsv(path_dist_tsv) %>%
  arrange(group.a, group.b) %>%
  select(group.a, group.b, euclidean.distance) %>%
  pivot_wider(names_from = group.b, values_from = euclidean.distance) %>%
  column_to_rownames(var = "group.a")

jaccard_hc <- as.dist(jaccard_dist_df) %>% hclust()
euclidean_hc <- as.dist(euclidean_dist_df) %>% hclust()

plot(
  jaccard_hc,
  
  # Label at same height
  hang = -1,
  main = 'primate14.chr6',
  xlab = 'Haplotype',
  ylab = 'Jaccard distance',
  sub = '',
  cex = 1.2
)


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
