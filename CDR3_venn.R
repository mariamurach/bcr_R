library(tidyverse)
#suppressMessages(library(ggpubr))
#library(ggpubr)
#library(rstatix)
library(tidyr)
#library(ggprism)
library(VennDiagram)
library(gridExtra)
 
source("~/universal.R")

samples <- read_csv("../samples.csv")
rownames(samples) <- samples$library
bcr.heavy <- read_csv("bcr_heavy.csv.gz")

plotVennCDR3 <- function(cell1,
                         cell1_cond,
                         cell2,
                         cell2_cond, main_font = 1, colA, colB)
{
  a <- bcr.heavy %>% filter(cell_type == cell1 & condition == cell1_cond)
  #b1a_WT_CDR3 <- b1a_WT %>% pull(CDR3aa) %>% unique
  a.seq <- a  %>% pull(CDR3aa) %>% unique
  
  b <- bcr.heavy %>% filter(cell_type == cell2 & condition == cell2_cond)
  #b1a_KO_CDR3 <- b1a_KO %>% pull(CDR3aa) %>% unique
  b.seq <- b %>% pull(CDR3aa) %>% unique
  ab <- intersect(a.seq, b.seq) %>% length
  a.perc <- round((100 - ab/length(a.seq) * 100), 1)
  b.perc <- round((100 - ab/length(b.seq) * 100), 1)
  a.cat <- paste0(a.perc, "%")
  b.cat <- paste0(b.perc, "%")
  pl <- venn.diagram(
    x = list(a.seq, b.seq),
    category.names = c(a.cat , b.cat),
    #filename = paste0('plots/', cell1, "_", cell1_cond, "_", cell2, "_", cell2_cond, '.png'),
    filename = NULL, 
    output=TRUE,
    print.mode=c("raw"),
    # # Output features
    imagetype="png" ,
    height = 480,
    width = 480,
    resolution = 300,
    compression = "lzw",
    # 
    # # Circles
    lwd = 2,
    lty = 'blank',
    fill = c(colA, colB),
    alpha = 0.5,
    # 
    # # Numbers
    cex =main_font,
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
    #cat.fontfamily = "ARIAL",
    cat.cex=main_font * 0.8,
    disable.logging = T,
    main.just = c(0.5, 0.5), 
    sub.just = c(0.5, 0.5), 
    sub.pos = c(180, 180),
    
    #sub.pos = c(0.5,0.5),
    ext.dist = -0.16,
    #ext.line.lwd = 0,
    #ext.dist = -0.27,
    ext.length = 0.83
  )
}

pl1 <- plotVennCDR3("B-1a", "WT", "B-1a", "KO", colA = vivid_colors[1], colB = vivid_colors[2])
pl2 <- plotVennCDR3("B-1a", "WT", "B-1b", "WT", colA = vivid_colors[1], colB = vivid_colors[3])
pl3 <- plotVennCDR3("B-1a", "WT", "B-1b", "KO", colA = vivid_colors[1], colB = vivid_colors[4])
pl4 <- plotVennCDR3("B-1a", "KO", "B-1b", "WT", colA = vivid_colors[2], colB = vivid_colors[3])
pl5 <- plotVennCDR3("B-1a", "KO", "B-1b", "KO", colA = vivid_colors[2], colB = vivid_colors[4])
pl6 <- plotVennCDR3("B-1b", "WT", "B-1b", "KO", colA = vivid_colors[3], colB = vivid_colors[4])
plots <- list()
plots[[6]] <- pl1
plots[[5]] <- pl2
plots[[3]] <- pl3
plots[[4]] <- pl4
plots[[2]] <- pl5
plots[[1]] <- pl6
m <- matrix(NA, 3, 3)
m[lower.tri(m, diag = T)] <- 1:6
all_no_legend <- grid.arrange(grobs = plots, layout_matrix = m)
ggsave(plot = all_no_legend, filename = "venns.toedit.pdf", w = 10, h = 10, unit = "in")

