library(ggplot2)
library(gggenes)
library(RColorBrewer)
library(ggpubr)

# molecule gene start end strand direction (reverse/-: start>end)

setwd("ab/spacers_new/gggenes/")

genes <-read.table("all.gggenes", header = TRUE, sep = "\t") # en el bak está el anterior de viruses

mis_colores <- c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
                 "#8C564B", "#C49C94", "#7F7F7F", "#C7C7C7", "#E377C2", "#F7B6D2", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5",
                 "#414451", "#FFFFFF")
genes$molecule <- factor(genes$molecule, levels = c("Prophage (2217)","CRISPR_I-Fa (225)","Choline_transporter (2671)",
                            "Trehalose_biosynthesis (1810)","Endonuclease (306)","Exonuclease (2071)",
                            "Restriction-Methylation_1 (41)","Restriction-Methylation_2 (129)","Restriction-Methylation_3 (324)",
                            "Nucleoid-associated_protein (148)","Asparagine_synthase (68)"))
gg <- ggplot(genes, aes(xmin = start, xmax = end, y = molecule, label = gene, fill = description2)) +
  geom_gene_arrow(show.legend = T, arrow_body_height = grid::unit(5, "mm"), 
                  arrowhead_width = grid::unit(7, "mm"), arrowhead_height = grid::unit(7, "mm")) +
  geom_gene_label(align = "centre", min.size = 2, height = grid::unit(4, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_manual(values = mis_colores) + ylab("Genomic island") +
  theme_genes() + theme (legend.title = element_blank(), legend.position = "none", legend.text = element_text(size = 7.2)) +
  guides(fill = guide_legend(nrow = 2))
print(gg)
pos <- ggplot() + theme_void() + geom_text(aes(0, 0, label = "Nucleotide position (bp)"))
pos <- cowplot::plot_grid(NULL, pos, nrow = 1, rel_widths = c(0.2, 1))

leg <- cowplot::get_legend(gg + theme(legend.position = "top"))
ggg <- cowplot::plot_grid(leg, gg, pos, ncol = 1, rel_heights = c(0.1, 1, 0.03))

setwd("Figures/def/")
pdf("fig2.pdf", width = 12, height = 7)
print(ggg)
dev.off()
