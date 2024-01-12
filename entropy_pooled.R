CDR3aa <- pool_CDR3_aa(bcr.heavy.v, "B-1b")

# Sample BCR or CDR3 sequences for two lists
list1_sequences <- CDR3aa %>% filter(condition == "WT") %>% pull(CDR3aa)
list2_sequences <- CDR3aa %>% filter(condition == "KO") %>% pull(CDR3aa)
# Function to calculate Shannon entropy
calculate_entropy <- function(sequence) {
  probabilities <- table(strsplit(sequence, NULL)) / nchar(sequence)
  entropy_value <- entropy::entropy(probabilities)
  return(entropy_value)
}
# Calculate entropy for each list
entropy_list1 <- sapply(list1_sequences, calculate_entropy)
entropy_list2 <- sapply(list2_sequences, calculate_entropy)
entropy_list_wt <- tibble(condition = "WT", entropies = entropy_list1, CDR3aa = names(entropy_list1)) %>%
  left_join(CDR3aa, by = c("condition" = "condition", "CDR3aa" = "CDR3aa"))
entropy_list_ko <- tibble(condition = "KO", entropies = entropy_list1, CDR3aa = names(entropy_list1)) %>%
  left_join(CDR3aa, by = c("condition" = "condition", "CDR3aa" = "CDR3aa"))


entropy_all <- bind_rows(entropy_list_wt, entropy_list_ko)
entropy_all$entropy_times_freq <- entropy_all$entropies*entropy_all$f_CDR3
entropy_all <- entropy_all %>% aggregate(entropy_times_freq ~ condition, sum)
entropy_all$condition <- factor(entropy_all$condition, levels = c("WT", "KO"))
entropy_all %>% ggplot(aes(x = condition, y = entropy_times_freq, fill = condition)) + 
  geom_bar(position = "dodge", stat = "identity") + theme_prism() + 
  scale_fill_manual(values = colo) + ggtitle("B-1b") + ylab("Shannon Entropy")
ggsave("plots/B-1b_entropy_pooled.png", w = 3.3, h = 2.5)