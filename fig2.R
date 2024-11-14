library(ggplot2)
library(gggenes)
library(plasmapR)
library(RColorBrewer)
#library(paletteer)
library(ggpubr)
library(stringr)
library(reshape2)
library(tidyverse)
library(dplyr)
library(cowplot)

# molecule gene start end strand direction (reverse/-: start>end)

na_strings <- ""
genes <-read.csv("all_freqs.gggenes", header = TRUE, sep = "\t", na = na_strings) # en el bak está el anterior de viruses

#getPalette <- colorRampPalette(brewer.pal(12, "Set3"))
#mis_colores <- getPalette(21)
#mis_colores[13] <- mis_colores[12]
#mis_colores[12] <- "#FFFFFF"
mis_colores <- c("#FFFFFF", "#414451", "#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", "#98DF8A", "#D62728", "#FF9896", 
                 "#9467BD", "#C5B0D5", "#8C564B", "#C49C94", "#7F7F7F", "#C7C7C7", "#E377C2", "#F7B6D2", "#BCBD22", "#DBDB8D",
                 "#17BECF", "#9EDAE5")

mis_colores <- c("#FF7F0E", "#FFBB78", "#1F77B4", "#AEC7E8", "#2CA02C", "#98DF8A", "#E377C2", "#F7B6D2", "#9467BD", "#C5B0D5",
                 "#BCBD22", "#DBDB8D")

mis_colores <- c("#FFB66C", "#71B260", "#b3b3a9", "#1F77B4", "#ec7374", "#E377C2", "#facde1", "#8158cb", "#8C564B", "#BCBD22")

genes <- genes %>% mutate(across('molecule', str_replace, "Prophage", "DgiS1 prophage"))
genes$molecule <- factor(genes$molecule, levels = c("Choline_transporter (2671)","DgiS1 prophage (2217)","Gao_Qat (2071)",
                            "PD-Lambda-5 (1810)","Restriction-Modification_1 (324)","Restriction-Modification_2 (306)","CRISPR_I-Fa (225)",
                            "Nucleoid-associated_protein (148)","Restriction-Modification_3 (129)",
                            "Asparagine_synthase (68)","Restriction-Modification_Septu (41)"))
gg <- ggplot(genes, aes(xmin = start, xmax = end, y = molecule, label = gene, fill = description)) + #, fill = gene
  geom_gene_arrow(show.legend = T, arrow_body_height = grid::unit(5, "mm"), 
                  arrowhead_width = grid::unit(6, "mm"), arrowhead_height = grid::unit(7, "mm")) +
  geom_gene_label(align = "centre", min.size = 6, height = grid::unit(5, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_manual(values = mis_colores, na.value = "white") + 
  ylab("Genomic island") + theme_genes() + 
  theme (legend.position = "none", legend.text = element_text(size = 12)) +
  guides(fill = guide_legend(nrow = 1)) + labs(fill = "Defense systems") 

pos <- ggplot() + theme_void() + geom_text(aes(0, 0, label = "Nucleotide position (bp)"))
pos <- cowplot::plot_grid(NULL, pos, nrow = 1, rel_widths = c(0.2, 1))

leg <- get_plot_component(gg + 
                            scale_fill_manual(values = mis_colores, na.translate = F) + 
                            theme(legend.position = "top") + 
                            labs(fill = "Defense systems"),
                            "guide-box-top")
ggg <- cowplot::plot_grid(leg, gg, pos, ncol = 1, rel_heights = c(0.1, 1, 0.03))
print(ggg)

# Heatmap
pa <- as.data.frame(readxl::read_excel("defense_islands.xlsx", sheet = "Heatmap"))
df <- data.frame(reshape2::melt(pa))
df[is.na(df)] <- 0

df2 <- as.data.frame(readxl::read_excel("defense_islands.xlsx", sheet = "Heatmap_2"))

df2$Islands <- factor(df2$Islands, levels = rev(pa$Islands))
#hm <- ggplot(df, aes(x = variable, y = Islands, fill = factor(value))) +
hm <- ggplot(df2, aes(x = variable, y = Islands, fill = value)) +
  geom_point(size = 6, shape = 21, color = "transparent") +
#  scale_fill_manual(values = c("transparent", "#7D1E6D")) +
  scale_fill_gradientn(colors = c("#E9BFDD", "#7D1E6D", "#25123C"), name = "% genomes", na.value = "transparent") +
  theme_minimal() +
  scale_x_discrete(position = "top") +
  xlab("Defense system") +
  ylab("Genomic island") +
  theme(axis.text.x = element_text(angle = 30, vjust = 0, hjust = 0),
        plot.margin = margin(r = 1, unit = "cm")
#        legend.position = "none"
        )
print(hm)

# Figure
fig <- ggarrange(ggg, NULL, hm, 
          labels = c("A", "", "B"), nrow = 3, heights = c(1, 0.035, 0.45))

pdf("fig2.pdf", width = 14, height = 10)
print(fig)
dev.off()

