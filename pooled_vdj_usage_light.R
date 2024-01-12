library(tidyverse)
suppressMessages(library(ggpubr))
library(ggpubr)
library(rstatix)
library(tidyr)
library(ggprism)

source("universal.R")

samples <- read_csv("samples.csv")
bcr<- read_csv("bcr.light.csv.gz")
bcr$condition <- factor(bcr$condition, levels = c("WT", "KO"))

all <- bcr
all$V <- gsub("\\-.*", "", all$V)
all$V <- gsub("\\*.*", "", all$V)

all$V <- gsub("S.*", "", all$V)
all$V <- gsub("IGKV", "VK", all$V)
all$V <- gsub("IGLV", "VL", all$V)
all$J <- gsub("\\*.*", "", all$J)
all$J <- gsub("IGLJ", "JL", all$J)
all$J <- gsub("IGKJ", "JK", all$J)
v_genes <- paste( "VK", seq(1:20),sep="")
v_genes <- c(v_genes, c("VL1", "VL2", "VL3"))
all$V <- factor(all$V, levels = v_genes)
get_vdj_plot <- function(cell.type, chain)
{
  genes <- all %>% filter(cell_type == cell.type, get(chain) != ".") %>%
    select(cell_type, condition, count, !!chain) %>% 
    group_by(cell_type, condition) %>%
    mutate(seq_n = sum(count)) %>%
    ungroup %>% 
    group_by(cell_type, condition, vdj_gene = get(chain))  %>% 
    dplyr::summarize(chain_count = sum(count), chain_freq = sum(count)/seq_n ) %>%
    ungroup %>% distinct %>% arrange(-chain_freq) %>% ungroup %>% drop_na()
  
  genes %>%  ggplot(aes(x = vdj_gene, y = chain_freq, fill = condition )) +
    geom_bar(stat = "identity", position = "dodge",  color = "black", linewidth = 1) +
    scale_fill_manual(values = vivid_colors) +
    scale_color_manual(values = vivid_colors) +
    theme_prism(base_size = 10) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
    ylab("% of all sequences") + xlab("")  +
    ggtitle(paste0(cell.type, ": ", chain, " genes")) +
    scale_y_continuous(labels = scales::percent)
}
get_vdj_plot("B-1a", "V")
ggsave("./plots/light_B-1a_V.pdf", h = 3, w = 9, dpi = 300)
ggsave("./plots/light_B-1a_V.png", h = 3, w = 9, dpi = 300)

get_vdj_plot("B-1a", "J")
ggsave("./plots/light_B-1a_J.png", h = 3, w = 3, dpi = 300)
ggsave("./plots/light_B-1a_J.pdf", h = 3, w = 3, dpi = 300)


get_vdj_plot("B-1b", "V")
ggsave("./plots/light_B-1b_V.pdf", h = 3, w = 9, dpi = 300)
ggsave("./plots/light_B-1b_V.png", h = 3, w = 9, dpi = 300)

get_vdj_plot("B-1b", "J")
ggsave("./plots/light_B-1b_J.png", h = 3, w = 3, dpi = 300)
ggsave("./plots/light_B-1b_J.pdf", h = 3, w = 3, dpi = 300)


