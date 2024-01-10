library(tidyverse)
suppressMessages(library(ggpubr))
library(ggpubr)
library(rstatix)
library(tidyr)
library(ggprism)

source("universal.R")

samples <- read_csv("samples.csv")
bcr.heavy <- read_csv("bcr_heavy.csv.gz")
bcr.heavy$condition <- factor(bcr.heavy$condition, levels = c("WT", "KO"))

all <- bcr.heavy
all$V <- gsub("\\-.*", "", bcr.heavy$V)
all$V <- gsub("S.*", "", all$V)
all$V <- gsub("IGHV", "VH", all$V)
all$D <- gsub("\\-.*", "", bcr.heavy$D)
all$D <- gsub("IGHD", "DH", all$D)
all$J <- gsub("\\*.*", "", bcr.heavy$J)
all$J <- gsub("IGHJ", "JH", all$J)
all$V <- factor(all$V, levels = c("VH1", "VH2", "VH3", "VH4", "VH5", "VH6", "VH7", "VH8", "VH9", "VH10", "VH11", "VH12", "VH13", "VH14", "VH15", "VH16"))


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
ggsave("./plots/B-1a_V.pdf", h = 3, w = 9, dpi = 300)
ggsave("./plots/B-1a_V.png", h = 3, w = 9, dpi = 300)

get_vdj_plot("B-1a", "J")
ggsave("./plots/B-1a_J.png", h = 3, w = 3, dpi = 300)
ggsave("./plots/B-1a_J.pdf", h = 3, w = 3, dpi = 300)

get_vdj_plot("B-1a", "D")
ggsave("./plots/B-1a_D.png", h = 3, w = 3, dpi = 300)
ggsave("./plots/B-1a_D.pdf", h = 3, w = 3, dpi = 300)


get_vdj_plot("B-1b", "V")
ggsave("./plots/B-1b_V.pdf", h = 3, w = 9, dpi = 300)
ggsave("./plots/B-1b_V.png", h = 3, w = 9, dpi = 300)

get_vdj_plot("B-1b", "J")
ggsave("./plots/B-1b_J.png", h = 3, w = 3, dpi = 300)
ggsave("./plots/B-1b_J.pdf", h = 3, w = 3, dpi = 300)

get_vdj_plot("B-1b", "D")
ggsave("./plots/B-1b_D.png", h = 3, w = 3, dpi = 300)
ggsave("./plots/B-1b_D.pdf", h = 3, w = 3, dpi = 300)

