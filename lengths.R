library(tidyverse)
library(ggpubr)
library(ggprism)
bcr <- read_csv("bcr_heavy.csv.gz")

bcr$cdr3_length <- nchar(bcr$CDR3nt)
bcr$condition <- factor(bcr$condition, levels = c("WT", "KO"))
source("universal.R")
CDR3aa <- pool_CDR3_nt(bcr, "B-1a")
CDR3aa$cdr3_length <- nchar(CDR3aa$CDR3nt)

CDR3aa <- CDR3aa %>% group_by(condition, cdr3_length) %>% dplyr::summarize(count = sum(c_CDR3), freq = sum(f_CDR3)) 
reps <- rep(CDR3aa$cdr3_length, CDR3aa$count)
cond <- rep(CDR3aa$condition, CDR3aa$count)
lengths <- tibble(cdr3_length = reps, condition = cond)
lengths %>% group_by(condition) %>% dplyr::summarize(mean(cdr3_length))

lengths %>% ggplot(aes(x = cdr3_length, fill = condition)) +
  geom_histogram(position = "dodge", bins = 50) +
  theme_prism() + facet_wrap(~condition) + ggtitle("B-1a") + scale_fill_manual(values = vivid_colors)+ xlab("CDR3_length")
ggsave("B1a_length.png", w = 6, h = 3.6)


CDR3aa <- pool_CDR3_nt(bcr, "B-1b")
CDR3aa$cdr3_length <- nchar(CDR3aa$CDR3nt)

CDR3aa <- CDR3aa %>% group_by(condition, cdr3_length) %>% dplyr::summarize(count = sum(c_CDR3), freq = sum(f_CDR3)) 
reps <- rep(CDR3aa$cdr3_length, CDR3aa$count)
cond <- rep(CDR3aa$condition, CDR3aa$count)
lengths <- tibble(cdr3_length = reps, condition = cond)

lengths %>% ggplot(aes(x = cdr3_length, fill = condition)) +
  geom_histogram(position = "dodge", bins = 50) +
  theme_prism() + facet_wrap(~condition) + ggtitle("B-1b")+ 
scale_fill_manual(values = vivid_colors) + xlab("CDR3_length")
ggsave("B1b_length.png", w = 6, h = 3.6)
lengths %>% group_by(condition) %>% dplyr::summarize(mean(cdr3_length))

