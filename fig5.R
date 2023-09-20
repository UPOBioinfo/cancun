library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(ggpubr)
library(ggbreak)
require(ggseqlogo)

setwd("/home/ajperez/Nextcloud/ncbidatasets/ab/spacers_new")

nsp <- read.csv("nspacers3.tsv", sep = "\t", header = T)
nsp$Array <- factor(nsp$Array, levels = c("ifb", "ifa", "cho", "tnp", "vir"))
nsp$N <- factor(nsp$N)
colores <- c("#1F78B4", "#FF7F00", "#E31A1C", "#6A3D9A", "#33A02C", "#FDBF6F", "#FB9A99", "#CAB2D6", "#B2DF8A")
rename(count(nsp, Array, N), Freq = n)
p1 <- ggplot(nsp, aes(x = Array, y = Frequency, fill = N)) +
  geom_boxplot() + stat_compare_means(aes(group = N), label = "p.format") +
  scale_fill_manual(values = c("#FF7F00", "#FDBF6F", "#1F78B4")) +
  ylab("No. of spacers per array") +
  scale_x_discrete(labels= c("I-Fb", "I-Fa (ssrA)", "I-Fa (cho)", "I-Fa (tnp)", "I-Fa (vir)")) +
  xlab("CRISPR array type") + labs(fill = "CRISPR\narray") +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 12), legend.position = "right",
        legend.text = element_text(size = 8), legend.title = element_text(size = 10))
print (p1)

##

nsp2 <- read.csv("nspacers_types.tsv", sep = "\t", header = T)
p2 <- ggplot(nsp2, aes(reorder(Type, -Frequency), Frequency, fill = reorder(Type, -Frequency))) +
  geom_bar(stat="identity",  fill = "#FB8072") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 12), legend.position = "none",
        legend.text = element_text(size = 8), legend.title = element_text(size = 10)) +
  ylab("Number of different spacers") +
  xlab("CRISPR array type(s)") +
  scale_y_continuous(trans = "log10")
  #scale_y_break(c(700, 2500), scales = 0.4, ticklabels = c(2500, 2600, 2700, 2800, 2900))
print(p2)

# Repeats
setwd("/home/ajperez/Nextcloud/ncbidatasets/ab/spacers_new/pam")

fig <- list()
cs1 = make_col_scheme(chars=c('A', 'T', 'C', 'G'), groups=c('gr1', 'gr2', 'gr3', 'gr4'), 
                      cols=c('#66C2A5', '#FC8D62', '#8DA0CB', '#E5C494'))

i <- 0
nombres <- c("I-Fa1", "I-Fa2", "I-Fb", "Vir1", "Vir2", "Cho1", "Cho2", "Tnp1", "Tnp2")
for (type in c("ifa1", "ifa2", "ifb", "vir1", "vir2", "bet1", "bet2", "tra1", "tra2")) {
  i = i + 1
  seqs <- readLines(paste0("rp_", type, ".seq"))
  fig[[type]] <- ggseqlogo(seqs, namespace = 'ATCG', col_scheme=cs1) + #method = 'p'
    ggtitle(nombres[i]) +
    theme(plot.margin = unit(c(-0.6, 1.5, -0.2, 0.6), "cm"), legend.position = "none", plot.title = element_text(hjust = 1.1, vjust = -10),
          text = element_text(size = 12), axis.text.x=element_blank())
  print(fig[[type]])
}

# add x axis
fig[["vir2"]] <- ggseqlogo(seqs, namespace = 'ATCG', col_scheme=cs1) + #method = 'p'
  ggtitle(nombres[i]) +
  theme(plot.margin = unit(c(-0.6, 1.5, 0, 0.6), "cm"), legend.position = "none", plot.title = element_text(hjust = 1.1, vjust = -10),
        text = element_text(size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Nucleotide position")

# Add annotates
fig[["ifb"]] <- fig[["ifb"]] +
  annotate('rect', xmin = 4.4, xmax = 5.55, ymin = -0.1, ymax = 2.1, alpha = 0.2, col='black', fill='blue') +
  annotate('rect', xmin = 6.4, xmax = 9.55, ymin = -0.1, ymax = 2.1, alpha = 0.2, col='black', fill='blue') +
  annotate('rect', xmin = 15.4, xmax = 19.55, ymin = -0.1, ymax = 2.1, alpha = 0.2, col='black', fill='blue')
fig[["ifa1"]] <- fig[["ifa1"]] +
  annotate('rect', xmin = 3.4, xmax = 4.55, ymin = -0.1, ymax = 2.1, alpha = 0.2, col='black', fill='yellow')
fig[["ifa2"]] <- fig[["ifa2"]] +
  annotate('rect', xmin = 3.4, xmax = 4.55, ymin = -0.1, ymax = 2.1, alpha = 0.2, col='black', fill='yellow')

p3 <- ggarrange(fig[["ifb"]], fig[["ifa1"]], fig[["bet1"]], 
          fig[["tra1"]], fig[["vir1"]], fig[["ifa2"]], fig[["bet2"]], fig[["tra2"]], fig[["vir2"]], 
          ncol = 1, heights = c(1,1,1,1,1,1,1,1,1.3))
print(p3)

# Output
fig1 <- ggarrange(p1, print(p2), 
                 labels = c("A", "B"), ncol = 1, heights = c(0.8, 1))
fig2 <- ggarrange(p3, labels = "C")
fig <- ggarrange (print(fig1), print(fig2), ncol = 2, width = c(0.5, 1))
print(fig)

setwd("/home/ajperez/Documentos/Articulos/ABA Project/ABA2 - Cancun/Figures/def/")
pdf("fig5.pdf", width = 12, height = 7)
print(fig)
dev.off()

