# R code for Steele et al., 2024 BioRxiv: https://doi.org/10.1101/2024.08.09.607352
# Compiled August 26, 2024

### Fisher's Exact Test
# Input file has five columns: category, subset count, subset total, genome count, genome total
Pl <- read.csv("counts.csv", row.names = 1, header= TRUE)
fisher.pval <- apply(Pl, 1, function(x) fisher.test(matrix(x,nr=2))$p.value)
fisher.pval.adj <- p.adjust(fisher.pval, method = "fdr")
View(fisher.pval)
View(fisher.pval.adj)
write.csv(fisher.pval.adj, "pvals.csv")


### Construct barplots
library(ggplot2)

COG_fractions <- read.csv("Pseudomonas_COG_fractions.csv")

GOG_barplot <- ggplot(COG_fractions, aes(Type, Percent, fill = COG)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = 60)) +
  scale_fill_manual(values = c("#015100","#0E700C","#248D21","#3FAA3C","#8DE58A",
                               "#010A60","#101971","#222B81", "#323C97","#3944AC","#525DC2","#6A75D8","#7E89F1","#98A1FC",
                               "#67001F","#840707","#B22424","#D94949","#EF6F65","#F78276","#FE9F95","#FEBCB5",
                               "#838383","#2D2C2C","#ACACAC")) +
  facet_wrap(~Strain) + theme_classic()
 GOG_barplot

social_fractions <- read.csv("Pseudomonas_social_fractions.csv")

social_barplot <- ggplot(social_fractions, aes(Type, Percent, fill = Category2)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = 60)) +
  scale_fill_manual(values = c("#2E2E2E","#4F4F4F","#757474","#FFFFFF","#969696","#B9B8B8","#2E2E2E","#DEDDDD")) +
  facet_grid(Category1~Strain) + theme_classic()
social_barplot 


### Construct heatmap
library(RColorBrewer)
library(reshape2)
library(ggplot2)

pident_table <- read.csv("Pseudomonas_predation_resistance_homologs.csv", row.names=1, na.strings="0")
data <- as.matrix(pident_table)
data[is.na(data)] <- 0
data_melt <- melt(data)

ggp <- ggplot(data_melt, aes(Var2, Var1, fill = value)) +                         
  geom_tile(color="black") + coord_fixed() + 
  theme(axis.text.x = element_text(angle = 60)) +
  scale_fill_gradient(low = "white", high = "black") + # #330066
  xlab("Genes important for predation resistance from Tn screen") + ylab("Strain") + labs(fill = "% identity")
ggp


# Construct tanglegrams
library(ape)
library(dendextend)
library(Quartet)

# Import Newick format tree files
# Prot_tree is a maximum likelihood phylogeny generated from an alignment for a single protein of interest
# Spec_tree is a maximum likelihood phylogeny from a concatenated alignment of 78 conserved proteins
prot_tree <- "prot.phy.treefile"
spec_tree <- "spec.phy.treefile"
prot_unrooted <- read.tree(prot_tree)
spec_unrooted <- read.tree(spec_tree)

prot_unrooted$tip.label <- as.character(prot_unrooted$tip.label)
spec_unrooted$tip.label <- as.character(spec_unrooted$tip.label)

# Calculate tree distances
quartet_distance <- QuartetPoints(spec_unrooted, prot_unrooted)
print(quartet_distance)
euclidean_distance <- treedist(spec_unrooted, prot_unrooted)
print(euclidean_distance)

# Don't specify outgroups, it's more trouble than it's worth

# Convert phylogeny to binary tree, then convert the binary tree to a dendrogram
prot_binary <- ape::multi2di(prot_unrooted)
prot_dist <- cophenetic(prot_binary)
prot_dist[is.na(prot_dist)] <- 0
prot_dist[is.infinite(prot_dist)] <- max(prot_dist[is.finite(prot_dist)])
prot_hclust <- hclust(as.dist(prot_dist))
prot_dend <- as.dendrogram(prot_hclust)

spec_binary <- ape::multi2di(spec_unrooted)
spec_dist <- cophenetic(spec_binary)
spec_dist[is.na(spec_dist)] <- 0
spec_dist[is.infinite(spec_dist)] <- max(spec_dist[is.finite(spec_dist)])
spec_hclust <- hclust(as.dist(spec_dist))
spec_dend <- as.dendrogram(spec_hclust)

# Create a dendlist
dl1 <- dendlist(spec_dend, prot_dend)

# Untangle the dendrograms to reduce crossing lines
dl1_untangled <- dl1 %>% untangle(method = "step2side") # method = c("labels", "ladderize", "random", "step1side", "step2side", "DendSer")

# Plot the untangled tanglegram
tanglegram(dl1_untangled, 
           sort = TRUE, 
           common_subtrees_color_lines = FALSE, 
           highlight_distinct_edges = FALSE, 
           highlight_branches_lwd = FALSE, 
           common_subtrees_color_branches = FALSE)


