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
bcr <- bcr %>% filter(CDR3aa != "out_of_frame", !str_detect(CDR3aa, "_"), !str_detect(CDR3aa, "\\?"))

bcr_b1a <- pool_CDR3_aa(bcr, "B-1a")
bcr_b1a$hydrophobicity <- hydrophobicity(bcr_b1a$CDR3aa)

bcr_summary_mean_b1a <- bcr_b1a %>% group_by(condition) %>% dplyr::summarize(hydrophobicity = mean(hydrophobicity))

bcr_summary_mean_b1a   %>%
  ggplot(aes(x = condition, y = hydrophobicity, color = condition)) + 
  geom_point(size = 4) + theme_prism() + 
  scale_color_manual(values = colo)  + ylab("Mean Hydrophobicity")+ ggtitle("B-1a")
ggsave("plots/mean_hydrophobicity_b1a.png", w = 5, h = 2.7)
ggsave("plots/mean_hydrophobicity_b1a.pdf", w = 5, h = 2.7)

bcr_summary_median_b1a <- bcr_b1a %>% group_by(condition) %>% dplyr::summarize(hydrophobicity = median(hydrophobicity))

bcr_summary_median_b1a  %>%
  ggplot(aes(x = condition, y = hydrophobicity, color = condition)) + 
  geom_point(size = 4) + theme_prism() + 
  scale_color_manual(values = colo)  + ylab("Median Hydrophobicity") + ggtitle("B-1a")
ggsave("plots/median_hydrophobicity_b1a.png", w = 5, h = 2.7)
ggsave("plots/median_hydrophobicity_b1a.pdf", w = 5, h = 2.7)

bcr_b1a  %>%
  ggplot(aes(x = condition, y = hydrophobicity, fill = condition)) + 
  geom_boxplot() + theme_prism() + 
  scale_fill_manual(values = colo) + ggtitle("B-1a") + 
  ylab("hydrophobicity")+ stat_compare_means() + coord_cartesian(ylim = c(-2.5, 2.5))
ggsave("plots/boxplot_hydrophobicity_b-1a.png", w = 3.5, h = 4)
ggsave("plots/boxplot_hydrophobicity_b-1a.pdf", w = 3.5, h = 4)


bcr_b1b <- pool_CDR3_aa(bcr, "B-1b")
bcr_b1b$hydrophobicity <- hydrophobicity(bcr_b1b$CDR3aa)

bcr_summary_mean_b1b <- bcr_b1b %>% group_by(condition) %>% dplyr::summarize(hydrophobicity = mean(hydrophobicity))

bcr_summary_mean_b1b   %>%
  ggplot(aes(x = condition, y = hydrophobicity, color = condition)) + 
  geom_point(size = 4) + theme_prism() + 
  scale_color_manual(values = colo)  + ylab("Mean Hydrophobicity") + ggtitle("B-1b")
ggsave("plots/mean_hydrophobicity_b1b.png", w = 5, h = 2.7)
ggsave("plots/mean_hydrophobicity_b1b.pdf", w = 5, h = 2.7)

bcr_summary_median_b1b <- bcr_b1b %>% group_by(condition) %>% dplyr::summarize(hydrophobicity = median(hydrophobicity))

bcr_summary_median_b1b  %>%
  ggplot(aes(x = condition, y = hydrophobicity, color = condition)) + 
  geom_point(size = 4) + theme_prism() + 
  scale_color_manual(values = colo)  + ylab("Median Hydrophobicity") + ggtitle("B-1b")
ggsave("plots/median_hydrophobicity_b1b.png", w = 5, h = 2.7)
ggsave("plots/median_hydrophobicity_b1b.pdf", w = 5, h = 2.7)

bcr_b1b  %>%
  ggplot(aes(x = condition, y = hydrophobicity, fill = condition)) + 
  geom_boxplot() + theme_prism() + 
  scale_fill_manual(values = colo) + ggtitle("B-1b") + 
  ylab("hydrophobicity")+ stat_compare_means(vjust = -1) + coord_cartesian(ylim = c(-2.5, 2.5))
ggsave("plots/boxplot_hydrophobicity_b-1b.png", w = 3.5, h = 4)
ggsave("plots/boxplot_hydrophobicity_b-1b.pdf", w = 3.5, h = 4)
