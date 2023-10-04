library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(ggpubr)
library(forcats)
library(ggVennDiagram)
library(readxl)

setwd("./")

inames <- as.data.frame(readxl::read_excel("islands.xlsx", sheet = "Islands"))
cancun <- read_excel("islands.xlsx", sheet = 3, col_names = TRUE, col_types = "text")

cancun <- cancun %>% filter(!if_all(p1:p6, is.na))
data <- data.frame(cancun %>% pivot_longer(p1:p6, names_to = "position", values_to = "gene") %>% 
                     group_by(position, gene) %>% 
                     dplyr::summarise(Frequency = n()))
data <- data %>% drop_na(gene)

colourCount <- nrow(unique(data %>% select(gene)))
getPalette <- colorRampPalette(c("black", brewer.pal(10, "Paired"), brewer.pal(9, "Set1")))
mis_colores <- getPalette(colourCount)
mis_colores[6] <- "darkorange"
mis_colores[10] <- "#E377C2"
mis_colores[12] <- "#EFC000"
mis_colores[14] <- mis_colores[5]
mis_colores[13] <- "#D04C4C"
mis_colores[5] <- "#00468B"
mi_breaks <- levels(factor(c("DgiS1","Choline_transporter","Trehalose_biosynthesis","CRISPR","Exonuclease","Endonuclease",
              "Restriction-methylation_1","Restriction-methylation_2","Restriction-methylation_3","Nucleoid-associated protein",
              "Asparagine_synthase","Other","Gene_core")))

# Hash
dmap <- hash::hash(inames$Gene, inames$Name)
data$gene <- sapply(strsplit(data$gene, ',\\s?'), function(x, na.rm = F) {
  paste(unlist(as.list(dmap)[x]), collapse = ',')
})
data$gene <- factor(data$gene, levels = inames$Name)

# Positions
coline = "#00008B"
siline = 1
tyline = 2
#sum(data %>% filter(position %in% "p2") %>% select(Frequency))
p1 <- ggplot(data, aes(x = position, y = Frequency, fill = gene)) +
  geom_bar(position="fill", stat="identity", width = 0.65) + 
  scale_fill_manual(values = mis_colores, breaks = levels(factor(data$gene))) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 0), text = element_text(size = 12), legend.position = "right",
        legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + 
  labs(fill = "Genomic island") + xlab("Position") +
  guides(colour = guide_legend(ncol1 = 2), fill = guide_legend(ncol = 2)) +
  scale_x_discrete(label = c("p1 (9,095)", "p2 (5,767)", "p3 (3,298)", "p4 (1,565)", "p5 (14)", "p6 (3)"))
print(p1)

alline <- 0.3
coline <- "red3"
d1 <- data.frame( x = c(1.34, 1.66, 1.66, 1.34), y = c(0.003, 0, 1, 1))
d2 <- data.frame( x = c(1.34+1, 1.66+1, 1.66+1, 1.34+1), y = c(0.5985, 0, 1, 1))
d3 <- data.frame( x = c(1.34+2, 1.66+2, 1.66+2, 1.34+2), y = c(0.5577, 0, 1, 1))
d4 <- data.frame( x = c(1.34+3, 1.66+3, 1.66+3, 1.34+3), y = c(0.9920, 0, 1, 1))
d5 <- data.frame( x = c(1.34+4, 1.66+4, 1.66+4, 1.34+4), y = c(0.8462, 0, 1, 1))
p1 <- p1 + 
  annotate("polygon", x = d1$x, y = d1$y, fill = coline, alpha = alline) + 
  annotate("polygon", x = d2$x, y = d2$y, fill = coline, alpha = alline) +
  annotate("polygon", x = d3$x, y = d3$y, fill = coline, alpha = alline) +
  annotate("polygon", x = d4$x, y = d4$y, fill = coline, alpha = alline) + 
  annotate("polygon", x = d5$x, y = d5$y, fill = coline, alpha = alline)

# Order
cancunf <- as.data.frame(readxl::read_excel("islands.xlsx", sheet = "Cancun(freq)"))

# Hash
cancunf$Order <- sapply(strsplit(cancunf$Order, ',\\s?'), function(x, na.rm = F) {
  paste(unlist(as.list(dmap)[x]), collapse = ',')
})

cancunf <- cancunf %>% filter(!Order == "None")
p2 <- ggplot(cancunf, aes(x = fct_rev(fct_inorder(Order)), y = Frequency, fill = factor(Color))) +
  geom_bar(stat="identity") + scale_fill_manual(values = c("darkorange", "black"), breaks = c("CRISPR-Cas I-Fa", "Prophage DgiS1")) +
  labs(fill = "Genomic island\nincluded") + xlab("Genomic islands ordered") +
  coord_flip() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 12), legend.position = "top",
        legend.text = element_text(size = 8), legend.title = element_text(size = 10))
print(p2)

# Venn diagram
setwd("ab/spacers_new")

p1virus <- readLines("p1virus.ab")
p1virus_i <- readLines("p1virus_incomplete.ab")
ifa <- readLines("../types/ifa.ab")
ifb <- readLines("../types/ifb.ab")
dvenn <- list(DgiS1 = p1virus, "DgiS1(incomplete)" = p1virus_i, "I-Fa" = ifa, "I-Fb" = ifb)

fvenn <- ggVennDiagram(dvenn, label_alpha = 0.5) + 
  scale_fill_gradient2(low = "white", high = "darkred") +
  scale_color_manual(values = c("black", "black", "black", "black")) +
  labs(fill = "No. genomes") +
  theme(legend.position = "top", legend.key.width = unit(1.1, "cm"))
print(fvenn)

# Output
fig2 <- ggarrange(p2, fvenn, 
          labels = c("B", "C"), nrow = 1, widths = c(1, 0.7))
fig <- ggarrange(p1, fig2, 
          labels = c("A", ""), nrow = 2)

setwd("Figures/def/")
pdf("fig3.pdf", width = 12, height = 10)
print(fig)
dev.off()

