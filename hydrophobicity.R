library(tidyverse)
library(ggpubr)
library(ggprism)
library(Peptides)
bcr <- read_csv("bcr_heavy.csv.gz")
colo <- c('#4363d8', '#f58231','#e6194b', '#3cb44b', '#ffe119',  '#911eb4',
          '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
          '#9a6324', '#fffac8', '#800000', '#aaffc3',
          '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
bcr$condition <- factor(bcr$condition, levels = c("WT", "KO"))
bcr$hydrophobicity <- hydrophobicity(bcr$CDR3aa)


bcr_summary_mean <- bcr %>% group_by(cell_type, condition) %>% dplyr::summarize(hydrophobicity = mean(hydrophobicity))

bcr_summary_mean  %>% filter(cell_type != "B2") %>%
  ggplot(aes(x = condition, y = hydrophobicity, color = condition)) + 
  geom_point(size = 4) + theme_prism() + 
  scale_color_manual(values = colo)  + ylab("Mean Hydrophobicity") + facet_wrap(~cell_type)
ggsave("plots/mean_hydrophobicity.png", w = 5, h = 2.7)

bcr_summary <- bcr %>% group_by(cell_type, condition) %>% dplyr::summarize(hydrophobicity = median(hydrophobicity))

bcr_summary  %>% filter(cell_type != "B2") %>%
  ggplot(aes(x = condition, y = hydrophobicity, color = condition)) + 
  geom_point(size = 4) + theme_prism() + 
  scale_color_manual(values = colo)  + ylab("Median Hydrophobicity") + facet_wrap(~cell_type)
ggsave("plots/median_hydrophobicity.png", w = 5, h = 2.7)


bcr  %>% filter(cell_type != "B2") %>%
  ggplot(aes(x = condition, y = hydrophobicity, fill = condition)) + 
  geom_boxplot() + theme_prism() + 
  scale_fill_manual(values = colo) + ggtitle("B-1a") + ylab("hydrophobicity") + facet_wrap(~cell_type) + stat_compare_means()
