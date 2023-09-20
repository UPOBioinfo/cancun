library(ggplot2)
library(gggenes)
library(RColorBrewer)

# molecule gene start end strand direction (reverse/-: start>end)

setwd("/home/ajperez/Nextcloud/ncbidatasets/ab/spacers_new/gggenes/")

genes <-read.table("all.gggenes", header = TRUE, sep = "\t")

#getPalette <- colorRampPalette(c("grey", "#fb8072"))
getPalette <- colorRampPalette(brewer.pal(12, "Set3"))
mis_colores <- getPalette(12)
mis_colores[13] <- mis_colores[12]
mis_colores[12] <- "#FFFFFF"
genes$molecule <- factor(genes$molecule, levels = c("P1virus (2217)","CRISPR_I-Fa (225)","Choline_transporter (2671)",
                            "Trehalose_biosynthesis (1810)","Endonuclease (306)","Exonuclease (2071)",
                            "Restriction-Methylation_1 (41)","Restriction-Methylation_2 (129)","Restriction-Methylation_3 (324)",
                            "Nucleoid-associated_protein (148)","Asparagine_synthase (68)"))
gg <- ggplot(genes, aes(xmin = start, xmax = end, y = molecule, label = gene, fill = unknown)) + #, fill = gene
  geom_gene_arrow(show.legend = T, arrow_body_height = grid::unit(5, "mm"), 
                  arrowhead_width = grid::unit(7, "mm"), arrowhead_height = grid::unit(7, "mm")) +
  geom_gene_label(align = "centre", min.size = 2, height = grid::unit(4, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  #scale_fill_manual(palette = "Set3") +
  scale_fill_manual(values = mis_colores) + ylab("Genomics island") +
  theme_genes() + theme (legend.title = element_blank(), legend.position = "top", legend.text = element_text(size = 7.2)) +
  guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))
print(gg)

setwd("/home/ajperez/Documentos/Articulos/ABA Project/ABA2 - Cancun/Figures/def/")
pdf("fig2.pdf", width = 12, height = 7)
print(gg)
dev.off()
