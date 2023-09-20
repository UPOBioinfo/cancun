require(ggplot2)
require(ggseqlogo)
library(gggenes)
require(ggpubr)
library(RColorBrewer)
library(forcats)
library(ggbreak)
library(ggrepel)
library(dplyr)

setwd("/home/ajperez/Nextcloud/ncbidatasets/ab/spacers_new/pam")

nombres <- c("I-Fa1", "I-Fa2", "I-Fb")
cs1 = make_col_scheme(chars=c('A', 'T', 'C', 'G'), groups=c('gr1', 'gr2', 'gr3', 'gr4'), cols=c('#66C2A5', '#FC8D62', '#8DA0CB', '#E5C494'))

x1 <- 1.46
x2 <- 1.435

# P1virus
#########
i <- 0
figl1 <- list()
for (type in c("ifa1", "ifa2", "ifb")) { # faltan los vir, xq están vacios
  i = i + 1
  sepp = 1
  if(type == "ifa2") { sepp = 0 } # align ifa2
  seqs <- readLines(paste0("pam_p1virus_", type, "_12.seq"))
  n <- length(seqs)
  figl1[[type]] <- ggseqlogo(seqs, namespace = 'ATCG', col_scheme=cs1) +
    theme(plot.margin = unit(c(-0.5, 2.1, 0, 0.6), "cm"), legend.position = "none", plot.title = element_text(hjust = x1 - 0.05, vjust = -11),
          text = element_text(size = 12), axis.text.x=element_blank()) + ggtitle(paste0(nombres[i], strrep(" ", sepp), "\n(", n, ")")) + 
    scale_x_continuous(breaks = 1:10, labels = -10:-1)
}

# Add X axis to I-Fb
figl1[[type]] <- ggseqlogo(seqs, namespace = 'ATCG', col_scheme=cs1) + #method = 'p'
  ggtitle(nombres[i]) +
  theme(plot.margin = unit(c(-0.5, 2.1, 0, 0.6), "cm"), legend.position = "none", plot.title = element_text(hjust = x2 - 0.05, vjust = -11),
        text = element_text(size = 12)) + ggtitle(paste0(nombres[i], " \n (", n, ")")) + 
  scale_x_continuous(breaks = 1:10, labels = -10:-1) +
  xlab("Nucleotide position")
ff1 <- ggarrange(figl1[["ifa1"]], figl1[["ifa2"]], figl1[["ifb"]], ncol = 1, hjust = -2, heights = c(1, 1, 1.05))
print(ff1)

# Add tittle
text <- "P1virus"
tgrob <- text_grob(text,size = 16)
titulo1 <- as_ggplot(tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))
ff1 <- ggarrange(titulo1, ff1, ncol = 1, heights = c(0.1, 1))

## Phigaros
###########
i <- 0
figl2 <- list()
for (type in c("ifa1", "ifa2", "ifb")) { # faltan los vir, xq están vacios
  i = i + 1
  seqs <- readLines(paste0("pam_phigaros_", type, "_12.seq"))
  n <- length(seqs)
  figl2[[type]] <- ggseqlogo(seqs, namespace = 'ATCG', col_scheme=cs1) +
    theme(plot.margin = unit(c(-0.5, 2.1, 0, 0.6), "cm"), legend.position = "none", plot.title = element_text(hjust = x1, vjust = -11),
          text = element_text(size = 12), axis.text.x=element_blank()) + ggtitle(paste0(nombres[i], "  \n(", n, ")")) + 
    scale_x_continuous(breaks = 1:10, labels = -10:-1)
}

# Add X axis to I-Fb
figl2[[type]] <- ggseqlogo(seqs, namespace = 'ATCG', col_scheme=cs1) + #method = 'p'
  ggtitle(nombres[i]) +
  theme(plot.margin = unit(c(-0.5, 2.1, 0, 0.6), "cm"), legend.position = "none", plot.title = element_text(hjust = x2, vjust = -11),
        text = element_text(size = 12)) + ggtitle(paste0(nombres[i], "   \n(", n, ")")) + 
  scale_x_continuous(breaks = 1:10, labels = -10:-1) +
  xlab("Nucleotide position")
ff2 <- ggarrange(figl2[["ifa1"]], figl2[["ifa2"]], figl2[["ifb"]], ncol = 1, hjust = -2, heights = c(1, 1, 1.05))
print(ff2)

# Add tittle
text <- "Remaining phages"
tgrob <- text_grob(text,size = 16)
titulo2 <- as_ggplot(tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))
ff2 <- ggarrange(titulo2, ff2, ncol = 1, heights = c(0.1, 1))

# Accumulations through p1virus
###############################
setwd("/home/ajperez/Nextcloud/ncbidatasets/ab/spacers_new/pam/")

#df <- read.csv("spacers_vs_ab00002_80_80_2.tsv", sep = "\t")
df <- read.csv("spacers_vs_p1virusffn_ditribution.tsv", sep = "\t")

labels <- df |> filter(Position %in% c(80815,66014,41670,25499,80678,6862,30980,40245,21899,79013,70145,46836,4191,5693,8742,10002,66668))
p <- ggplot(df, aes(x = Position, y = Frequency)) + #, label = paste0(Label, "\n", N)
  geom_bar(stat="identity", width = 5, color = "#C36B74") + # scale_x_continuous(breaks=seq(1, 81980, 10000)) +
  #geom_text(aes(label = paste0(Label, "\n", N)), lineheight = 0.75) + # vjust = 0, hjust = 0
  geom_label_repel(aes(label = paste0(Label, "\n(", N, ")")), data = labels, 
                   size = 4,
                   force_pull = 0, 
                   nudge_x = 0.5,
                   box.padding = 0.5,
                   nudge_y = 0.5,
                   min.segment.length = 0, # draw all lines no matter how short
                   segment.size = 0.2,
                   #segment.curvature = -0.1,
                   segment.ncp = 3,
                   segment.angle = 45,
                   label.size = NA,
                   fill = NA) +   ylab("No. of protospacers per gene") +
  theme_classic() +  
  theme (axis.line.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) # x4=izq(-: atrás), x2=der(-:alante)
#print(p)

# add gene browser
genes <-read.table("../gggenes/p1virus.gggenes", header = TRUE, sep = "\t")

g <- ggplot(genes, aes(xmin = start, xmax = end, y = molecule, label = gene, fill = factor(Spacers))) + #, fill = unknown)) +
  geom_gene_arrow(show.legend = F, arrow_body_height = grid::unit(5, "mm"), arrowhead_width = grid::unit(7, "mm"), arrowhead_height = grid::unit(7, "mm")) +
  geom_gene_label(align = "centre", min.size = 2, height = grid::unit(4, "mm")) +
  scale_fill_manual(values =  c("#FFFFFF", "#C36B74"))+
  theme_genes() + theme (legend.position = "none", axis.title.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"),
                         axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank())
p <- p + coord_cartesian(xlim = c(0, 81981))
g <- g + coord_cartesian(xlim = c(0, 81981))
pg <- list(p, g)
fig2 <- ggarrange(plotlist = pg, ncol = 1, align = "v", heights = c(1, 0.2))
print(fig2)

# Freq PAM
##########
setwd("/home/ajperez/Nextcloud/ncbidatasets/ab/spacers_new/pam/")
cc <- read.csv("pam_freq_p1virus_cds.tsv", header = F, sep = "\t")
colnames(cc) <- c("Frequency", "PAM")
mis_colores <- c(brewer.pal(12, "Set3"))
mis_colores[2] <- "#FFED6F"
pam <- ggplot (cc, aes(x = fct_rev(fct_reorder(PAM, Frequency)), y = Frequency, fill = PAM)) +
  geom_col() + theme_minimal() +
  #scale_y_break(c(2000, 10000)) +
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = mis_colores) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 12), legend.position = "none",
        legend.text = element_text(size = 8), legend.title = element_text(size = 10)) +
  xlab("PAM sequence")
pam
# Output
#########
fig1 <- ggarrange(ff1, print(pam), ff2, nrow = 1, labels = c("A", "B", "C"))
fig <- ggarrange(fig1, NULL, fig2, labels = c("", "D", ""), ncol = 1, heights = c(1, 0.1, 1))
print(fig)

setwd("/home/ajperez/Documentos/Articulos/ABA Project/ABA2 - Cancun/Figures/def/")
pdf("fig7.pdf", width = 12, height = 8)
print(fig)
dev.off()

