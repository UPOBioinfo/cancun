library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(ggpubr)
library(ggbreak)
require(ggseqlogo)
library(stringr)
library(ggplotify)
library(UpSetR)
library(ComplexUpset)

nsp <- read.csv("nspacers3.tsv", sep = "\t", header = T)
nsp$Array <- factor(nsp$Array, levels = c("ifb", "ifa", "cho", "tnp", "vir"))
nsp$N <- factor(nsp$N)
colores <- c("#1F78B4", "#FF7F00", "#E31A1C", "#6A3D9A", "#33A02C", "#FDBF6F", "#FB9A99", "#CAB2D6", "#B2DF8A")
rename(count(nsp, Array, N), Freq = n)

#kruskal.test(Frequency ~ Array, data = nsp_stat)
#pairwise.wilcox.test(nsp$Frequency, nsp$Array, p.adjust.method = "BH")

nsp$N <- str_replace(nsp$N, "1", "C1a")
nsp$N <- str_replace(nsp$N, "2", "C2a")

p1 <- ggplot(nsp, aes(x = Array, y = Frequency, fill = N)) +
  geom_boxplot() + stat_compare_means(aes(group = N), label = "p.format") +
  scale_fill_manual(values = c("#FF7F00", "#FDBF6F", "#1F78B4")) +
  ylab("No. of spacers per array") +
  scale_x_discrete(labels= c("I-Fb", expression(paste("I-Fa (", italic("ssrA"), ")")), "I-Fa (choline)", "I-Fa (other)", "I-Fa (PPTOP)")) +
  xlab("CRISPR array type") + labs(fill = "CRISPR\narray") +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 12), legend.position = "right",
        legend.text = element_text(size = 8), legend.title = element_text(size = 10)) +
  geom_segment(x = 1.8, xend = 2.2, y = 170, yend = 170, size = 1) +
  geom_segment(x = 2.8, xend = 3.2, y = 170, yend = 170, size = 1) +
  geom_segment(x = 3.8, xend = 4.2, y = 170, yend = 170, size = 1) +
  geom_segment(x = 4.8, xend = 5.2, y = 170, yend = 170, size = 1) +
  geom_text(label = "***", x = 2, y = 172, size = 6) +
  geom_text(label = "ns", x = 3, y = 178, size = 4) +
  geom_text(label = "**", x = 4, y = 172, size = 6) +
  geom_text(label = "**", x = 5, y = 172, size = 6)
print (p1)

##
# New panel upset
lista <- list()
lista[["I-Fb"]] <- readLines("ifb.sp")
lista[["C1a (ssrA)"]] <- readLines("ifa1.sp")
lista[["C2a (ssrA)"]] <- readLines("ifa2.sp")
lista[["C1a (choline)"]] <- readLines("bet1.sp")
lista[["C2a (choline)"]] <- readLines("bet2.sp")
lista[["C1a (other)"]] <- readLines("tra1.sp")
lista[["C2a (other)"]] <- readLines("tra2.sp")
lista[["C1a (PPTOP)"]] <- readLines("vir1.sp")
lista[["C2a (PPTOP)"]] <- readLines("vir2.sp")

p2 <- upset(fromList(lista), name = "Combinations",
            rev(c("I-Fb", "C1a (ssrA)", "C2a (ssrA)", 
                  "C1a (choline)", "C2a (choline)", 
                  "C1a (other)", "C2a (other)", "C1a (PPTOP)", "C2a (PPTOP)")),
            queries=list(
              upset_query(set='I-Fb', fill='#1F78B4'),
              upset_query(set='C1a (ssrA)', fill='#FF7F00'),
              upset_query(set='C2a (ssrA)', fill='#FDBF6F'),
              upset_query(set='C1a (choline)', fill='#FF7F00'),
              upset_query(set='C2a (choline)', fill='#FDBF6F'),
              upset_query(set='C1a (other)', fill='#FF7F00'),
              upset_query(set='C2a (other)', fill='#FDBF6F'),
              upset_query(set='C1a (PPTOP)', fill='#FF7F00'),
              upset_query(set='C2a (PPTOP)', fill='#FDBF6F'),
              upset_query(intersect = c('I-Fb', 'C1a (ssrA)'), color='#D83842', fill='#D83842', only_components=c('Intersection size')),
              upset_query(intersect = c('C1a (ssrA)', 'C2a (ssrA)'), color='#D83842', fill='#D83842', only_components=c('Intersection size')),
              upset_query(intersect = c('I-Fb', 'C2a (ssrA)'), color='#D83842', fill='#D83842', only_components=c('Intersection size')),
              upset_query(intersect = c('I-Fb', 'C1a (ssrA)', 'C2a (ssrA)'), color='#D83842', fill='#D83842', only_components=c('Intersection size')),
              upset_query(intersect = c('I-Fb', 'C1a (other)'), color='#D83842', fill='#D83842', only_components=c('Intersection size')),
              upset_query(intersect = c('C1a (ssrA)', 'C2a (ssrA)', 'C1a (choline)'), color='#D83842', fill='#D83842', only_components=c('Intersection size')),
              upset_query(intersect = c('I-Fb', 'C1a (ssrA)', 'C1a (choline)'), color='#D83842', fill='#D83842', only_components=c('Intersection size')),
              upset_query(intersect = c('C1a (ssrA)', 'C2a (PPTOP)'), color='#D83842', fill='#D83842', only_components=c('Intersection size')),
              upset_query(intersect = c('I-Fb', 'C1a (ssrA)', 'C1a (other)'), color='#D83842', fill='#D83842', only_components=c('Intersection size')),
              upset_query(intersect = c('C2a (ssrA)', 'C1a (choline)'), color='#D83842', fill='#D83842', only_components=c('Intersection size'))
            ),
            base_annotations = list(
              'Intersection size' = (
                intersection_size(
                  text=list(
                    vjust = 0.1,
                    hjust = -0.1,
                    angle = 30,
                    size = 2
                  ),
                  bar_number_threshold = 1,  # show all numbers on top of bars
                  width = 0.85,   # reduce width of the bars
                )
                # add some space on the top of the bars
                + scale_y_continuous(expand=expansion(mult=c(0, 0.05)), limits = c(0, 3000))
                + theme(
                  text = element_text(size = 7.5),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  axis.line = element_line(colour='black')
                ) + ylab('Number of spacers in each combination')
              )
            ),
            stripes = upset_stripes(
              geom = geom_segment(linewidth = 12),  # make the stripes larger
              colors=c('grey95', 'white')
            ),
            # to prevent connectors from getting the colorured
            # use `fill` instead of `color`, together with `shape='circle filled'`
            matrix = intersection_matrix(
              geom = geom_point(
                shape='circle filled',
                size = 2.5,
                stroke = 0.45
              )
            ),
            set_sizes=(
              upset_set_size(geom=geom_bar(width = 0.7))
              + theme(
                text = element_text(size = 7.5),
                axis.line.x = element_line(colour = 'black'),
                axis.ticks.x = element_line()
              ) + ylab('Number of total spacers')
            ),
            sort_sets = F,
            sort_intersections = 'descending',
            themes=upset_modify_themes(
              list('intersections_matrix'=theme(axis.text.y=element_text(size = 7)))
            )
)
print(p2)

# Repeats
fig <- list()
cs1 = make_col_scheme(chars=c('A', 'T', 'C', 'G'), groups=c('gr1', 'gr2', 'gr3', 'gr4'), 
                      cols=c('#66C2A5', '#FC8D62', '#8DA0CB', '#E5C494'))

i <- 0
nombres <- c("C1a \n(ssrA)", "C2a \n(ssrA)", "I-Fb", "C1a \n(choline)", "C2a \n(choline)", 
             "C1a \n(other)", "C2a \n(other)", "C1a \n(PPTOP)", "C2a \n(PPTOP)")
for (type in c("ifa1", "ifa2", "ifb", "bet1", "bet2", "tra1", "tra2", "p-p1", "p-p2")) {
  i = i + 1
  seqs <- readLines(paste0("rp_", type, ".seq"))
  fig[[type]] <- ggseqlogo(seqs, namespace = 'ATCG', col_scheme=cs1) + #method = 'p'
    ggtitle(nombres[i]) +
    theme(plot.margin = unit(c(-0.6, 1.5, -0.2, 0.6), "cm"), legend.position = "none", plot.title = element_text(hjust = 1.1, vjust = -11),
          text = element_text(size = 8), axis.text.x=element_blank())
  print(fig[[type]])
}

# add x axis
fig[["p-p2"]] <- ggseqlogo(seqs, namespace = 'ATCG', col_scheme = cs1) + #method = 'p'
  ggtitle(nombres[i]) +
  theme(plot.margin = unit(c(-0.6, 1.5, 0, 0.6), "cm"), legend.position = "none", plot.title = element_text(hjust = 1.1, vjust = -11),
        text = element_text(size = 8), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
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
          fig[["tra1"]], fig[["p-p1"]], fig[["ifa2"]], fig[["bet2"]], fig[["tra2"]], fig[["p-p2"]], 
          ncol = 1, heights = c(1,1,1,1,1,1,1,1,1.3))
print(p3)

# Output
fig1 <- ggarrange(print(p1), NULL, print(p2), 
                 labels = c("A", "", "B"), ncol = 1, heights = c(1, 0.01, 1))
fig2 <- ggarrange(p3, labels = "C")
fig <- ggarrange (print(fig1), print(fig2), ncol = 2, width = c(0.5, 1))
print(fig)

pdf("fig5.pdf", width = 12, height = 7)
print(fig)
dev.off()

