library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(ggpubr)
library(forcats)
library(ggVennDiagram)

setwd("/home/ajperez/Documentos/Articulos/ABA Project/ABA2 - Cancun/")

inames <- as.data.frame(readxl::read_excel("islands.xlsx", sheet = "Islands"))
cancun <- as.data.frame(readxl::read_excel("islands.xlsx", sheet = "Cancun"))

cancun <- cancun %>% filter(!if_all(p1:p5, is.na))
data <- data.frame(cancun %>% pivot_longer(p1:p5, names_to = "position", values_to = "gene") %>% 
                     group_by(position, gene) %>% 
                     dplyr::summarise(Frequency = n()))
colourCount <- nrow(unique(data %>% select(gene)))
getPalette <- colorRampPalette(c("black", brewer.pal(10, "Paired"), brewer.pal(9, "Set1")))
mis_colores <- getPalette(colourCount)
mis_colores[6] <- "darkorange"
mis_colores[10] <- "#E377C2"
mis_colores[12] <- "#EFC000"
mis_colores[14] <- mis_colores[5]
mis_colores[13] <- mis_colores[14]
mis_colores[5] <- "#00468B"
mi_breaks <- levels(factor(c("P1virus","Choline_transporter","Trehalose_biosynthesis","CRISPR","Exonuclease","Endonuclease",
              "Restriction-methylation_1","Restriction-methylation_2","Restriction-methylation_3","Nucleoid-associated protein",
              "Asparagine_synthase","Other","Gene_core")))

# Hash
dmap <- hash::hash(inames$Gene, inames$Name)
data$gene <- sapply(strsplit(data$gene, ',\\s?'), function(x, na.rm = F) {
  paste(unlist(as.list(dmap)[x]), collapse = ',')
})
data$gene <- factor(data$gene, levels = inames$Name)

# Positions
p1 <- ggplot(data, aes(x = position, y = Frequency, fill = gene)) +
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = mis_colores, breaks = levels(factor(data$gene))) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 12), legend.position = "right",
        legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + 
  labs(fill = "Genomic\nisland") + xlab("Position") +
  guides(colour = guide_legend(ncol1 = 2), fill = guide_legend(ncol = 2))
print(p1)

# Order
cancunf <- as.data.frame(readxl::read_excel("islands.xlsx", sheet = "Cancun(freq)"))

# Hash
cancunf$Order <- sapply(strsplit(cancunf$Order, ',\\s?'), function(x, na.rm = F) {
  paste(unlist(as.list(dmap)[x]), collapse = ',')
})

p2 <- ggplot(cancunf, aes(x = fct_rev(fct_inorder(Order)), y = Frequency, fill = factor(Color))) +
  geom_bar(stat="identity") + scale_fill_manual(values = c("darkorange", "black"), breaks = c("CRISPR", "P1virus")) +
  labs(fill = "Genomic island\nincluded") + xlab("Genomic islands ordered") +
  coord_flip() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 12), legend.position = "top",
        legend.text = element_text(size = 8), legend.title = element_text(size = 10))
print(p2)

# Venn diagram
setwd("/home/ajperez/nc/ncbidatasets/ab/spacers_new")

p1virus <- readLines("p1virus.ab")
p1virus_i <- readLines("p1virus_incomplete.ab")
ifa <- readLines("../types/ifa.ab")
ifb <- readLines("../types/ifb.ab")
dvenn <- list(P1virus = p1virus, "P1virus(incomplete)" = p1virus_i, "I-Fa" = ifa, "I-Fb" = ifb)

fvenn <- ggVennDiagram(dvenn, label_alpha = 0.5) + 
  scale_fill_gradient2(low = "white", high = "darkred") +
  scale_color_manual(values = c("black", "black", "black", "black")) +
  labs(fill = "No. genomes") +
#  geom_rect(aes(xmin = 0, xmax = 0.4, ymin = 1, ymax = 58), fill = "white") +
  theme(legend.position = "top", legend.key.width = unit(1.1, "cm"))
print(fvenn)

# Output
fig2 <- ggarrange(p2, fvenn, 
          labels = c("B", "C"), nrow = 1, widths = c(1, 0.7))
fig <- ggarrange(p1, fig2, 
          labels = c("A", ""), nrow = 2)

#fig <- ggarrange(p1, p2, 
#          labels = c("A", "B"), nrow = 1, widths = c(0.9, 1))

setwd("/home/ajperez/Documentos/Articulos/ABA Project/ABA2 - Cancun/Figures/def/")
pdf("fig3.pdf", width = 12, height = 10)
print(fig)
dev.off()

