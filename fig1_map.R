#devtools::install_github("bradyajohnston/plasmapr")

#install.packages("remotes")
#remotes::install_github("BradyAJohnston/plasmapR", force = T)
#install.packages("stringi")

devtools::install("plasmapR")

library(tidyverse)
library(plasmapR)

dat <- read_tsv("pptop.gff3", comment = "#", col_names = FALSE) %>% 
  select(type = 3, start = 4, end = 5, direction = 7, attribute = 9) %>% 
  filter(!type %in% "gene") %>% 
  mutate(direction = ifelse(direction == "+", 1, -1), 
         type = ifelse(type == "region", "source", type),
         name = ifelse(type == "source",
                       str_match(attribute, "Name=(\\w+);")[,2],
                       str_match(attribute, "locus_tag=(\\w+);")[,2])) %>% 
  rownames_to_column("index")

plot_plasmid(dat, dat$name[1]) + theme(legend.position = "right")

dat[dat$type == "CDS", ] |> 
  plot_plasmid(name = "PPTOP")

dat |> 
  plot_plasmid(name = "pETM-20")

