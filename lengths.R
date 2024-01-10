library(tidyverse)
library(ggpubr)
library(ggprism)
bcr <- read_csv("bcr_heavy.csv.gz")

bcr$cdr3_length <- nchar(bcr$CDR3nt)
bcr$condition <- factor(bcr$condition, levels = c("WT", "KO"))

colo <- c('#4363d8', '#f58231','#e6194b', '#3cb44b', '#ffe119',  '#911eb4',
          '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
          '#9a6324', '#fffac8', '#800000', '#aaffc3',
          '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
bcr %>% filter(cell_type == "B-1a") %>% ggplot(aes(x = cdr3_length, fill = condition)) +
  geom_histogram(position = "dodge", bins = 50) +
  theme_prism() + facet_wrap(~condition) + ggtitle("B-1a") + scale_fill_manual(values = colo)
ggsave("B1a_length.png", w = 6, h = 3.6)
bcr %>% filter(cell_type == "B-1b") %>% ggplot(aes(x = cdr3_length, fill = condition)) +
  geom_histogram(position = "dodge", bins = 50) +
  theme_prism() + facet_wrap(~condition) + ggtitle("B-1b")+ scale_fill_manual(values = colo)
ggsave("B1b_length.png", w = 6, h = 3.6)
