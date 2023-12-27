library(tidyverse)
#suppressMessages(library(ggpubr))
#library(ggpubr)
#library(rstatix)
library(tidyr)
#library(ggprism)
library(VennDiagram)
 
source("~/universal.R")

samples <- read_csv("../samples.csv")
rownames(samples) <- samples$library
bcr.heavy <- read_csv("bcr_heavy.csv.gz")


a <- bcr.heavy %>% filter(cell_type == "B-1a" & condition == "WT")
#b1a_WT_CDR3 <- b1a_WT %>% pull(CDR3aa) %>% unique
a.seq <- a  %>% pull(CDR3aa) %>% unique

b <- bcr.heavy %>% filter(cell_type == "B-1a" & condition == "KO")
#b1a_KO_CDR3 <- b1a_KO %>% pull(CDR3aa) %>% unique
b.seq <- b %>% pull(CDR3aa) %>% unique
ab <- intersect(a.seq, b.seq) %>% length
a.perc <- round((100 - ab/length(a.seq) * 100), 1)
b.perc <- round((100 - ab/length(b.seq) * 100), 1)
a.cat <- paste(a.perc, "%")
b.cat <- paste(b.perc, "%")
pl <- venn.diagram(
  x = list(a.seq, b.seq),
  category.names = c(a.cat , b.cat),
  #filename = 'plots/B-1a_WT_KO.png',
  filename = NULL, 
  output=TRUE,
  print.mode=c("raw"),
  # # Output features
  imagetype="png" ,
  height = 480 ,
  width = 480 ,
  resolution = 300,
  compression = "lzw",
  # 
  # # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#f5b67a", "#7575fa"),
  # 
  # # Numbers
  cex = 1.1,
  fontface = "bold",
  fontfamily = "sans",
  # 
  # # Set names
  #cat.cex = 0.6,
  # cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(180, -180),
  # cat.dist = c(0.055, 0.055),
  # cat.fontfamily = "sans",
  # rotation = 1
  cat.fontfamily = "ARIAL",
  cat.cex=1.1,

)


grid::grid.draw(pl)
dev.off()
