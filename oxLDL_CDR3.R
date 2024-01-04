library(tidyverse)
suppressMessages(library(ggpubr))
library(ggpubr)
library(rstatix)
library(tidyr)
library(ggprism)

source("~/universal.R")
samples <- read_csv("samples.csv")
bcr.heavy <- read_csv("bcr_heavy.csv.gz")
bcr.heavy$condition <- factor(bcr.heavy$condition, levels = c("WT", "KO"))

CDR3aa <- bcr.heavy %>% 
  select(cell_type, condition,CDR3aa, count, C) %>% 
  group_by(cell_type, condition) %>%
  mutate(seq_n = sum(count)) %>%
  ungroup %>% 
  group_by(cell_type, condition, CDR3aa, C)  %>% 
  summarize(c_CDR3 = sum(count), f_CDR3 = sum(count)/seq_n ) %>%
  ungroup %>% distinct %>% arrange(-f_CDR3) %>% ungroup
CDR3_oxLDL <- CDR3aa %>% filter(CDR3aa == "CMRYGNYWYFDVW", cell_type != "B2", C == "IGHM")
###limit to IgM 
CDR3_oxLDL %>% ggplot(aes(x = condition, y = f_CDR3, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge",  color = "black", linewidth = 1) +
  scale_fill_manual(values = vivid_colors) +
  scale_color_manual(values = vivid_colors) +
  theme_prism(base_size = 10) %+replace% 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1), 
        title = element_text(size = 10),
        strip.text = element_text(
          size = 10, color = "dark green")) + 
  facet_wrap(~cell_type) +
  ylab("% of all sequences") + xlab("") + ggtitle("CDR3AA: CMRYGNYWYFDVW")
  