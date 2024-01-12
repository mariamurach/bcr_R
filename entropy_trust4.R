library(entropy)
library(tidyverse)
library(ggpubr)
library(ggprism)
library(bio3d)
bcr <- read_csv("bcr_heavy.csv.gz")
source("funcs.R")
source("universal.R")

cell = "B-1a"
bcr.heavy <- bcr
bcr.heavy$V <- gsub("\\-.*", "", bcr.heavy$V)
bcr.heavy$V <- gsub("S.*", "", bcr.heavy$V)
bcr.heavy$V <- gsub("IGHV", "VH", bcr.heavy$V)

bcr_b1a <- pool_CDR3_aa(dat = bcr, cell)
lin_b1a <- BuildBCRlineage(bcr_b1a)
clon_b1a <- getClonality(lin_b1a)
createDir("workspace")
save.image(paste0("workspace/" , cell, "_clonality.RData"))