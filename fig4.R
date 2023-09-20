library(dplyr)
library(ggplot2)
library(ggtree)
#library(ggbreak)
library(ggpubr)
library(stringr)
library(ape)
library(RColorBrewer)
library(tidyr)
library(grid)
library(ggflags)
library(Polychrome)

setwd("/home/ajperez/nc/ncbidatasets/ab/max")

mdata <- read.csv("metadata_ab_def_caudovirus.tsv", sep = "\t", header = T)
mfreq <- rename(count(mdata, MLST), Frequency = n)
mlst <- mfreq %>% filter(Frequency >= 50) %>% pull(MLST)

#mdata <- mdata %>% filter(MLST %in% mlst) %>% unite("CRISPR", Ambiguous, I.Fab, P1virus, sep = "/")
mdata <- mdata %>% filter(MLST %in% mlst) %>% unite("CRISPR", I.Fab, P1virus, sep = "/")
mdata <- rename(count(mdata, MLST, CRISPR), Frequency = n)
mdata$CRISPR <- str_replace(mdata$CRISPR, "//", "/")
mdata$CRISPR <- str_replace(mdata$CRISPR, "//", "/")
mdata$CRISPR <- str_replace(mdata$CRISPR, "^/", "")
mdata$CRISPR <- str_replace(mdata$CRISPR, "/$", "")
mdata$CRISPR <- str_replace(mdata$CRISPR, ",", "/")
mdata$MLST = paste("ST", mdata$MLST, sep='')
mdata$MLST <- str_replace(mdata$MLST, "ST-", "Unknown")
mdata$MLST <- factor(mdata$MLST, levels = c(paste0("ST", 0:10000), "Unknown"))

mi_breaks <- levels(factor(mdata$CRISPR))
mi_breaks[1] <- "No CRISPR"

colores <- c("azure4", "darkorange", "yellow3", "darkorchid", "darkblue", "#1F78B4", "lightblue", "black")

# Figure
p1 <- ggplot(mdata, aes(x = MLST, y = Frequency, fill = CRISPR)) +
  geom_bar(position="stack", stat="identity") + 
#  scale_y_break(c(500, 3500)) + scale_y_break(c(4000, 5500)) + 
  scale_fill_manual(values = colores, labels = mi_breaks) + #, breaks = levels(factor(mdata$CRISPR))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 12), legend.position = "none",
        legend.text = element_text(size = 8))
print(p1)
lp1 <- cowplot::get_legend(p1 + theme(legend.position = "right") + labs(fill = "Genomic islands"))

# Zoom
mdata2 <- mdata %>% filter(MLST != "ST2")
p2 <- ggplot(mdata2, aes(x = MLST, y = Frequency, fill = CRISPR)) +
  geom_bar(position="stack", stat="identity") + 
  #  scale_y_break(c(500, 3500)) + scale_y_break(c(4000, 5500)) + 
  scale_fill_manual(values = colores) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 12), legend.position = "none",
        legend.text = element_text(size = 8))
print(p2)

# Phylogeny
setwd("/home/ajperez/Nextcloud/ncbidatasets/ab/realphy")

# Phylo1
########
tree1 <- read.tree("tree.tree") # def
tree1 <- root(tree1, outgroup = "ab07408", resolve.root = TRUE) # distance / 10
heatmap_tabla <- read.table("metadata_st79.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = FALSE) # heatmap_def_core.tsv | clade_fig6_metadata.tsv | cancun.tsv;  row.names = 1,
heatmap_tabla2 <- read.table("metadata_st79.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = FALSE, row.names = 1)[1:2]
heatmap_tabla$Source <- factor(heatmap_tabla$Source, levels = c("J. Craig Venter Institute", "Oswaldo Cruz Institute"))

phylo1 <- ggtree(tree1, layout = 'rectangular', branch.length = 'variable', size = 1) %<+% heatmap_tabla + # none | variable
  geom_tiplab2(size = 0) 
#phylo1 <- ggtree(tree1, layout = 'rectangular', branch.length = 'variable', size = 1) + geom_tiplab2(size = 3, angle = 0) # see labels

phylo1 <- phylo1 +
  geom_tippoint(aes(color = Source, x = x + 0.00001), size = 3) +
  scale_color_manual(values = c("blue", "yellow3"), na.translate = F)

thispalette <- c("darkorange", "grey", "#377EB8", "#4DAF4A", "black")
phylo1 <- gheatmap(phylo1, heatmap_tabla2, offset = 0,  width = 0.3, font.size = 0, colnames_angle = 90, hjust = 1, colnames = F) +
  scale_fill_manual(values = thispalette, na.translate = F) +
  labs(fill = "Genome")
print(phylo1)

# Phylo2 (with flags)
#####################
tree <- read.tree("ppcrispr/data/CP087336/polymorphisms_move.phy_phyml_tree.txt") # def
tree <- root(tree, outgroup = "CP087336", resolve.root = TRUE)
heatmap_table <- read.table("metadata_ppcrispr.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = FALSE) # heatmap_def_core.tsv | clade_fig6_metadata.tsv | cancun.tsv;  row.names = 1,
heatmap_table2 <- read.table("metadata_ppcrispr.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = FALSE, row.names = 1)[c(3, 4, 5, 6, 7)]
#heatmap_table$Source <- factor(heatmap_table$Source, levels = c("Canada", "Thailand", "China", "Germany"))

phylo <- ggtree(tree, layout = 'rectangular', branch.length = 'variable', size = 1) %<+% heatmap_table + # none | variable
  geom_tiplab2(size = 0)

#phylo <- ggtree(tree, layout = 'rectangular', branch.length = 'variable', size = 1) + geom_tiplab2(size = 6, angle = 0) # see labels
phylo_hide <- ggtree(tree, layout = 'rectangular', branch.length = 'variable', size = 0)

#allScales1 <- c("Canada" = "#6F6F74", "Thailand" = "#BD0A36", "China" = "#F07E6E", "Germany" = "#000000", 
#               "Australia" = "#265DAB", "Colombia" = "#DECF3F", "Malaysia" = "#F17CB0", "Myanmar" = "#4BAF73", "Netherlands" = "#FCB11C",
#               "USA" = "#BC99C7", "Viet Nam" = "#9C755F")
allScales1 <- c("Canada" = "brown", "Thailand" = "cyan4", "China" = "red", "Germany" = "#000000", "Japan" = "white",
                "Malaysia" = "purple", "Myanmar" = "green2", "USA" = "#BC99C7")
allScales2 <- c("ST16" = "#FF7F0E", "ST32" = "#1F77B4", "ST150" = "#8C564B", "ST152" = "#D62728", "ST237" = "#9467BD", "ST626" = "#2CA028")

ph1 <- phylo +
  #geom_tippoint(aes(color = Source, x = x + 0.00005), shape = 18, size = 10, alpha = 0.9) +
  geom_tippoint(aes(fill = factor(MLST, levels = c("ST16", "ST32", "ST150", "ST152", "ST237", "ST626"))), shape = 23, size = 3.2) +
  #scale_color_manual(values = allScales1, name = "Source", na.translate = F) +
  scale_fill_manual(values = allScales2, name = "MLST", na.translate = F) +
  geom_flag(aes(country = Country, x = x + 0.002), na.rm = T, size = 5.5) +
  theme(legend.position = "none")
#print(ph1)
lph1 <- cowplot::get_legend(ph1 + theme(legend.position = "right") + labs(fill = "MLST"))

c <- c("ca","cn","de","mm","th");
p <- c("Canada","China","Germany","Myanmar","Thailand")
v <- 1:5
d <- data.frame(c, p, v)
paises <- ggplot(d, aes(x = p, y = v, country = c)) + geom_flag() + scale_country(labels = d$p, name = "Country") + theme(legend.position = "none") 
lpaises <- cowplot::get_legend(paises + theme_void() + theme(legend.position = "right"))

legends <- cowplot::plot_grid(lpaises, lph1, ncol = 1, rel_heights = c(0.1, 0.1))
#legends

phylo2 <- cowplot::plot_grid(ph1, legends, nrow = 1, rel_widths = c(0.7, 0.3))
print(phylo2)

# Phylo3
########
tree_3 <- read.tree("pp_completos/data/CP087336/polymorphisms_move.phy_phyml_tree.txt")
#tree <- root(tree, outgroup = "ab03783", resolve.root = TRUE)
heatmap_table_3 <- read.table("metadata_pp.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = FALSE) # heatmap_def_core.tsv | clade_fig6_metadata.tsv | cancun.tsv;  row.names = 1,
heatmap_table2_3 <- read.table("metadata_pp.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = FALSE, row.names = 1)[c(2, 3, 4, 5)]
#heatmap_table$Source <- factor(heatmap_table$Source, levels = c("Canada", "Thailand", "China", "Germany"))

phylo_3 <- ggtree(tree_3, layout = 'rectangular', branch.length = 'variable', size = 1) %<+% heatmap_table_3 + # none | variable
  geom_tiplab2(size = 0)
#phylo_3 <- ggtree(tree_3, layout = 'rectangular', branch.length = 'variable', size = 1) + geom_tiplab2(size = 6, angle = 0) # see labels
  
phylo_hide_3 <- ggtree(tree_3, layout = 'rectangular', branch.length = 'variable', size = 0)

mlsts <- readLines("metadata_pp.mlst")
c55 = createPalette(55,  c("#ff0000", "#00ff00", "#0000ff"))
c55 <- as.vector(c55)
c55[c(7, 9, 27, 28, 32, 40)] <- c("#FF7F0E","#1F77B4","#2CA02C","#D62728","#9467BD","#8C564B")
ph1_3 <- phylo_3 +
  #geom_tippoint(aes(color = Source, x = x + 0.00005), shape = 18, size = 10, alpha = 0.9) +
  geom_tippoint(aes(fill = factor(MLST, levels = mlsts)), shape = 23, size = 2) +
  #scale_color_manual(values = allScales1, name = "Source", na.translate = F) +
  scale_fill_manual(values = c55, name = "MLST", na.translate = F) +
  theme(legend.position = "none", legend.text = element_text(size = 6), legend.key.size = unit(0.5, 'cm'))
print(ph1_3)
lph1_3 <- cowplot::get_legend(ph1_3 + theme(legend.position = "right") + labs(fill = "MLST"))

n2_3 <- heatmap_table2_3 %>% select(P1virus)
ph2_3 <- gheatmap(phylo_hide_3, n2_3, offset = -10,  width = 450, font.size = 0, hjust = 1, colnames = F) +
  scale_fill_manual(values = c("#000000", "azure4"), name = "P1virus", na.translate = F) +
  theme(legend.position = "none")
#print(ph2_3)
lph2_3 <- cowplot::get_legend(ph2_3 + theme(legend.position = "right"))

n3_3 <- heatmap_table2_3 %>% select(Defense_system)
ph3_3 <- gheatmap(phylo_hide_3, n3_3, offset = -10,  width = 450, font.size = 0, colnames_angle = 90, hjust = 1, colnames = F) +
  scale_fill_manual(values = c("darkorange", "#FFE5CC", "#5CB85C", "#D54773"), name = "CRISPR-Cas", na.translate = F) +
  theme(legend.position = "none")
#print(ph3_3)
lph3_3 <- cowplot::get_legend(ph3_3 + theme(legend.position = "right"))

n4_3 <- heatmap_table2_3 %>% select(R.M)
ph4_3 <- gheatmap(phylo_hide_3, n4_3, offset = 0,  width = 450, font.size = 0, colnames_angle = 90, hjust = .5, colnames = F) +
  scale_fill_manual(values = c("#8074A8"), name = "R-M(p-p)", na.translate = F) +
  theme(legend.position = "none")
#print(ph4_3)
lph4_3 <- cowplot::get_legend(ph4_3 + theme(legend.position = "right"))

n5_3 <- heatmap_table2_3 %>% select(Coverage)
ph5_3 <- gheatmap(phylo_hide_3, n5_3, offset = -10,  width = 400, font.size = 0, colnames_angle = 90, hjust = 0, colnames = F) +
  scale_fill_continuous(low = "white", high = "darkred", name = "Coverage (%)") +
  geom_rect(aes(xmin = 0, xmax = 1.4, ymin = 1, ymax = 165), fill = "white") +
  theme(plot.margin = unit(c(5.5, 0, 5.5, 0), "pt"), legend.position = "none")
#print(ph5_3)
lph5_3 <- cowplot::get_legend(ph5_3 + theme(legend.position = "right"))

# Legend and final figure
n_3 <- heatmap_table2_3 %>% select(P1virus, Defense_system, R.M)
ph234_3 <- gheatmap(phylo_hide_3, n_3, offset = -1,  width = 300, font.size = 0, hjust = 1, colnames = F) +
  scale_fill_manual(values = c("darkorange", "#FFE5CC", "#5CB85C", "#D54773", "black", "azure4", "#8074A8"), name = "P1virus", na.translate = F) +
  geom_rect(aes(xmin = 0, xmax = 0.4, ymin = 1, ymax = 165), fill = "white") +
  theme(legend.position = "none")
print(ph234_3)

#legends1 <- cowplot::plot_grid(lpaises, lph1, ncol = 1, rel_heights = c(0.1, 0.1))
legends2_3 <- cowplot::plot_grid(lph2_3, lph3_3, lph4_3, lph5_3, ncol = 1, rel_heights = c(0.25, 0.35, 0.18, 0.45), align = "v")
legends_3 <- cowplot::plot_grid(legends2_3, lph1_3, ncol = 2)
#legends
phylo2_3 <- cowplot::plot_grid(ph1_3, ph234_3, ph5_3, NULL, legends_3, nrow = 1, rel_widths = c(0.4, 0.2, 0.08, 0.06, 0.58))
print(phylo2_3)

# Output
########
p1p2 <- ggarrange(p1, p2, labels = c("A", "B"), ncol = 1)
ab <- ggarrange(p1p2, lp1, phylo1, labels = c("", "", "C"), nrow = 1, widths = c(1, 0.5, 1))
cd <- ggarrange(phylo2_3, phylo2, labels = c("D", "E"), ncol = 2, widths = c(1, 0.5))
fig <- ggarrange(ab, cd, nrow = 2, heights = c(1, 0.75))

setwd("/home/ajperez/Documentos/Articulos/ABA Project/ABA2 - Cancun/Figures/def/")
pdf("fig4.pdf", width = 12, height = 10)
print(fig)
dev.off()
