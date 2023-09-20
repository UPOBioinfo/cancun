library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggpubr)

setwd("/home/ajperez/nc/ncbidatasets/ab/spacers_new")
t <- read.csv("variants_spacers2.tsv", sep = "\t") # variants_spacers2.tsv para ifa1 the same as tra1 vir1 bet1
t2 <- t %>% select("Phage_cluster", "fa_ifa1", "fa_ifa2", "fa_ifb", "fas_ifa1", "fas_ifa2", "fas_ifb", "fasp_ifa1", "fasp_ifa2", "fasp_ifb")
#colnames(t2) <- c("Phage", "I-Fa1 (no. variants)", "I-Fa2 (no. variants)", "I-Fb (no. variants)", "I-Fa1 (no. genomes)", "I-Fa2 (no. genomes)", "I-Fb (no. genomes)",
#                  "I-Fa1 (no. spacers)", "I-Fa2 (no. spacers)", "I-Fb (no. spacers)")

t <- t %>% select("Phage_cluster", "fr_ifa1", "fr_ifa2", "fr_ifb", "frs_ifa1", "frs_ifa2", "frs_ifb", "frsp_ifa1", "frsp_ifa2", "frsp_ifb")
#colnames(t) <- c("Phage", "I-Fa1 (% variants)", "I-Fa2 (% variants)", "I-Fb (% variants)", "I-Fa1 (% genomes)", "I-Fa2 (% genomes)", "I-Fb (% genomes)",
#                 "I-Fa1 (% spacers)", "I-Fa2 (% spacers)", "I-Fb (% spacers)")

t[,c(8,9,10)] <- t[,c(8,9,10)]*10
df <- melt(t, id.vars = "Phage_cluster", variable.name = "Type", value.name = "Count")
df2 <- melt(t2, id.vars = "Phage_cluster", variable.name = "Type", value.name = "Count")

#getPalette <- colorRampPalette(brewer.pal(9, "Set2")) # getPalette(9)
colores = c("#FF7F00", "#FDBF6F", "#1F78B4", "#FF7F00", "#FDBF6F", "#1F78B4", "#FF7F00", "#FDBF6F", "#1F78B4")
types1 <- c("fr_ifa1" = "I-Fa1", "fr_ifa2" = "I-Fa2", "fr_ifb" = "I-Fb",
            "frs_ifa1" = "I-Fa1", "frs_ifa2" = "I-Fa2", "frs_ifb" = "I-Fb",
            "frsp_ifa1" = "I-Fa1", "frsp_ifa2" = "I-Fa2", "frsp_ifb" = "I-Fb")
df <- df %>% mutate(Phage_cluster = factor(Phage_cluster), 
                    Phage_cluster = factor(Phage_cluster, levels = rev(levels(Phage_cluster))))
p1 <- ggplot(df, aes(Count, Phage_cluster, fill = Type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Type, nrow = 1, labeller = labeller(Type = types1)) +
  scale_fill_manual(values = colores) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_blank(), legend.position = "none", 
        plot.margin = unit(c(0.8, 0, 0, 0), "cm")) +
  ylab("Phages") + xlab("Percentage")
print(p1)

# Output
f1 <- ggarrange(p1)
fig1 <- f1 + annotate('segment', x = 0.1, xend = 0.39, y = 0.96, yend = 0.96, size = 1.2) +
  annotate('text', x = 0.24, y = 0.98, label = "% of recognized variants of the virus") +
  annotate('segment', x = 0.4, xend = 0.69, y = 0.96, yend = 0.96, size = 1.2) +
  annotate('text', x = 0.54, y = 0.98, label = "% of genomes with spacers against this phage") +
  annotate('segment', x = 0.70, xend = 0.99, y = 0.96, yend = 0.96, size = 1.2) +
  annotate('text', x = 0.84, y = 0.98, label = "% of different spacers")

setwd("/home/ajperez/Documentos/Articulos/ABA Project/ABA2 - Cancun/Figures/def/")
pdf("fig6.pdf", width = 12, height = 8)
print(fig1)
dev.off()

# Suppl. Fig.
types2 <- c("fa_ifa1" = "I-Fa1", "fa_ifa2" = "I-Fa2", "fa_ifb" = "I-Fb",
            "fas_ifa1" = "I-Fa1", "fas_ifa2" = "I-Fa2", "fas_ifb" = "I-Fb",
            "fasp_ifa1" = "I-Fa1", "fasp_ifa2" = "I-Fa2", "fasp_ifb" = "I-Fb")
df2 <- df2 %>% mutate(Phage_cluster = factor(Phage_cluster), 
                    Phage_cluster = factor(Phage_cluster, levels = rev(levels(Phage_cluster))))
p2 <- ggplot(df2, aes(Count, Phage_cluster, fill = Type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Type, nrow = 1, scales = "free_x", labeller = labeller(Type = types2)) +
  scale_fill_manual(values = colores) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1), legend.position = "none", 
        plot.margin = unit(c(0.8, 0, 0, 0), "cm")) +
   ylab("Phages")
print(p2)

f2 <- ggarrange(p2)
fig2 <- f2 + annotate('segment', x = 0.1, xend = 0.39, y = 0.96, yend = 0.96, size = 1.2) +
  annotate('text', x = 0.24, y = 0.98, label = "No. of recognized variants of the virus") +
  annotate('segment', x = 0.4, xend = 0.69, y = 0.96, yend = 0.96, size = 1.2) +
  annotate('text', x = 0.54, y = 0.98, label = "No. of genomes with spacers against this phage") +
  annotate('segment', x = 0.70, xend = 0.99, y = 0.96, yend = 0.96, size = 1.2) +
  annotate('text', x = 0.84, y = 0.98, label = "No. of different spacers")

setwd("/home/ajperez/Documentos/Articulos/ABA Project/ABA2 - Cancun/Figures/def/")
pdf("supplfigs3.pdf", width = 12, height = 8)
print(fig2)
dev.off()