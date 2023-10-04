library(ggplot2)
library(ggpubr)
library(pheatmap)
library(dplyr)
library(grid)

setwd("ab/cancun/phigaro")

virus <- readLines("../../spacers_new/virus_100_phage.id")
virus <- c("DgiS1", "DgiS1(incomplete)", "PPTOP", virus)
strains <- readLines("../../strains.ab")
gnmatrix <- read.csv("phage_matrix_p1virus.tsv", header = T, sep = '\t', row.names = 1)
gnmatrix <- gnmatrix %>% rename("DgiS1" = "P1virus")
gnmatrix <- gnmatrix %>% rename("DgiS1(incomplete)" = "P1virus.incomplete.")
gnmatrix <- gnmatrix %>% rename("PPTOP" = "Phage.plasmid")

gnmatrix <- gnmatrix %>% select(any_of(virus))
gnmatrix <- gnmatrix %>% filter(row.names(gnmatrix) %in% strains)
gnmatrix <- as.matrix(gnmatrix)

metadata <- read.csv("metadata_ab_phigaro.tsv", sep = "\t", header = T, row.names = 1)
metadata <- subset(metadata, rownames(metadata) %in% rownames(gnmatrix))
metadata <- metadata %>% select(CRISPR.Cas)
colnames(metadata) <- "CRISPR-Cas"

colores = list("CRISPR-Cas" = c("I-Fa" = "darkorange", "I-Fa,I-Fb" = "#E8DE9D", "I-Fb" = "#1F78B4", "ambiguous" = "#FFE5CC", "No_CRISPR" = "grey"))
p <- pheatmap(gnmatrix, show_rownames = 0, show_colnames = 1, treeheight_col = 0, legend = 0, angle_col = 90, 
              cluster_rows = T, cluster_cols = T, cex = 1, treeheight_row = 0, annotation_row = metadata,
              fontsize = 12, color = c("#F9F2EA", "#804AD5"), 
             annotation_colors = colores)
print(p)

setwd("Figures/def/")
pdf("supplfigs2.pdf", width = 12, height = 7)
grid.text("Genomes", x = 0.89, y = 0.55, rot = 90)
print(p)
dev.off()
