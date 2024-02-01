library(entropy)
library(tidyverse)
library(ggpubr)
library(ggprism)
library(bio3d)
bcr <- read_csv("bcr_heavy.csv.gz")
bcr.heavy.v <- bcr
bcr.heavy.v$V <- gsub("\\-.*", "", bcr.heavy.v$V)
bcr.heavy.v$V <- gsub("S.*", "", bcr.heavy.v$V)
bcr.heavy.v$V <- gsub("IGHV", "VH", bcr.heavy.v$V)

CDR3aa <- bcr.heavy.v %>% filter(cell_type == "B-1a") %>%
  select(cell_type, condition,CDR3aa, CDR3nt, count) %>% 
  group_by(cell_type, condition) %>%
  mutate(seq_n = sum(count)) %>%
  ungroup %>% 
  group_by(cell_type, condition, CDR3aa)  %>% 
  dplyr::summarize(c_CDR3 = sum(count), f_CDR3 = sum(count)/seq_n ) %>%
  ungroup %>% distinct %>% arrange(-f_CDR3)

CDR3aa.V <- bcr.heavy.v %>% filter(cell_type == "B-1a") %>%
  select(cell_type, condition,CDR3aa, CDR3nt, count, V) %>% 
  group_by(cell_type, condition) %>%
  mutate(seq_n = sum(count)) %>%
  ungroup %>% 
  group_by(cell_type, condition, CDR3aa, V)  %>% 
  dplyr::summarize(c_CDR3 = sum(count), f_CDR3 = sum(count)/seq_n ) %>%
  ungroup %>% distinct %>% arrange(-f_CDR3) %>% ungroup


CDR3_props <- CDR3aa.V %>% select(CDR3aa, condition, f_CDR3, c_CDR3) %>% group_by(condition) %>%  unique %>% ungroup
table(CDR3_props$condition, CDR3_props$f_CDR3 > 0.01)
seq_clones_WT <- sum((CDR3_props[CDR3_props$f_CDR3 > 0.01 & CDR3_props$condition == "WT",])$f_CDR3)
seq_clones_KO <- sum((CDR3_props[CDR3_props$f_CDR3 > 0.01 & CDR3_props$condition == "KO",])$f_CDR3)

# Sample BCR or CDR3 sequences for two lists
list1_sequences <- CDR3_props %>% filter(condition == "WT") %>% pull(CDR3aa)
list2_sequences <- CDR3_props %>% filter(condition == "KO") %>% pull(CDR3aa)

# Function to calculate Shannon entropy
calculate_entropy <- function(sequence) {
  probabilities <- table(strsplit(sequence, NULL)) / nchar(sequence)
  entropy_value <- entropy::entropy(probabilities)
  return(entropy_value)
}

colo <- c('#4363d8', '#f58231','#e6194b', '#3cb44b', '#ffe119',  '#911eb4',
          '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
          '#9a6324', '#fffac8', '#800000', '#aaffc3',
          '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
# Calculate entropy for each list
entropy_list1 <- sapply(list1_sequences, calculate_entropy)
entropy_list2 <- sapply(list2_sequences, calculate_entropy)
entropy_list_wt <- tibble(condition = "WT", entropies = entropy_list1)
entropy_list_ko <- tibble(condition = "KO", entropies = entropy_list2)
entropy_all <- bind_rows(entropy_list_wt, entropy_list_ko)
entropy_all <- entropy_all %>% aggregate(entropies ~ condition, sum)
entropy_all$condition <- factor(entropy_all$condition, levels = c("WT", "KO"))
entropy_all %>% ggplot(aes(x = condition, y = entropies, fill = condition)) + 
  geom_bar(position = "dodge", stat = "identity") + theme_prism() + 
  scale_fill_manual(values = colo) + ggtitle("B-1a") + ylab("Shannon Entropy")
ggsave("plots/B-1a_entropy.png", w = 3.3, h = 2.5)


CDR3aa <- bcr.heavy.v %>% filter(cell_type == "B-1b") %>%
  select(cell_type, condition,CDR3aa, CDR3nt, count) %>% 
  group_by(cell_type, condition) %>%
  mutate(seq_n = sum(count)) %>%
  ungroup %>% 
  group_by(cell_type, condition, CDR3aa)  %>% 
  dplyr::summarize(c_CDR3 = sum(count), f_CDR3 = sum(count)/seq_n ) %>%
  ungroup %>% distinct %>% arrange(-f_CDR3)

CDR3aa.V <- bcr.heavy.v %>% filter(cell_type == "B-1b") %>%
  select(cell_type, condition,CDR3aa, CDR3nt, count, V) %>% 
  group_by(cell_type, condition) %>%
  mutate(seq_n = sum(count)) %>%
  ungroup %>% 
  group_by(cell_type, condition, CDR3aa, V)  %>% 
  dplyr::summarize(c_CDR3 = sum(count), f_CDR3 = sum(count)/seq_n ) %>%
  ungroup %>% distinct %>% arrange(-f_CDR3) %>% ungroup


CDR3_props <- CDR3aa.V %>% select(CDR3aa, condition, f_CDR3, c_CDR3) %>% group_by(condition) %>%  unique %>% ungroup

# Sample BCR or CDR3 sequences for two lists
list1_sequences <- CDR3_props %>% filter(condition == "WT") %>% pull(CDR3aa)
list2_sequences <- CDR3_props %>% filter(condition == "KO") %>% pull(CDR3aa)

# Function to calculate Shannon entropy
calculate_entropy <- function(sequence) {
  probabilities <- table(strsplit(sequence, NULL)) / nchar(sequence)
  entropy_value <- entropy::entropy(probabilities)
  return(entropy_value)
}

colo <- c('#4363d8', '#f58231','#e6194b', '#3cb44b', '#ffe119',  '#911eb4',
          '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
          '#9a6324', '#fffac8', '#800000', '#aaffc3',
          '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
# Calculate entropy for each list
entropy_list1 <- sapply(list1_sequences, calculate_entropy)
entropy_list2 <- sapply(list2_sequences, calculate_entropy)
entropy_list_wt <- tibble(condition = "WT", entropies = entropy_list1)
entropy_list_ko <- tibble(condition = "KO", entropies = entropy_list2)
entropy_all <- bind_rows(entropy_list_wt, entropy_list_ko)
entropy_all <- entropy_all %>% aggregate(entropies ~ condition, sum)
entropy_all$condition <- factor(entropy_all$condition, levels = c("WT", "KO"))
entropy_all %>% ggplot(aes(x = condition, y = entropies, fill = condition)) + 
  geom_bar(position = "dodge", stat = "identity") + theme_prism() + 
  scale_fill_manual(values = colo) + ggtitle("B-1b") + ylab("Shannon Entropy")
ggsave("plots/B-1b_entropy.png",  w = 3.3, h = 2.5 )

