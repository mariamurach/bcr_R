library(tidyverse)
source("~/universal.R")
library(ggprism)
library(ggpubr)
bcr_light <- read_csv("bcr.light.csv.gz")
bcr_light$condition <- factor(bcr_light$condition, levels = c("WT", "KO"))

st.Ig <- bcr_light%>% 
  group_by(sample,condition, cell_type) %>%
  mutate(est_clonal_exp_norm = frequency/sum(frequency)) %>%  
  dplyr::filter(C != ".", cell_type != "B2") %>%
  group_by(sample,condition,cell_type, C) %>% 
  dplyr::summarise(Num.Ig = sum(est_clonal_exp_norm))  %>% ungroup

barplot_igs <- ggplot(st.Ig,aes(x = condition, y = Num.Ig, fill = C))+
  geom_bar(stat = "identity",position="fill",alpha = 0.8)+
  theme_bw()+
  labs(y = "Normalized Ig Abundance",fill = "Ig") + facet_wrap(~cell_type) +
  scale_fill_manual(values = c("#FD8008", "#0000FF")) +
  theme_prism()

st.Ig.long <- st.Ig %>% pivot_wider(names_from = C, values_from = Num.Ig)
st.Ig.long %>% group_by(cell_type, condition) %>% summarize(kappa_to_lambda = IGKC/IGLC)
ggsave("plots/lambda_kappa_ratio.png", barplot_igs, w = 5.3, h = 5.3, dpi = 300, unit = "in")


all_Bs <- st.Ig %>% mutate(C = factor(C)) %>% dplyr::filter(C != ".", cell_type != "B2") %>%
  ggplot(aes(x = condition, y = Num.Ig, fill = condition))+ 
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  labs(y = "Normalized Ig Abundance",fill = "Ig")+
  scale_fill_manual(values =  c("#0000FF", "#FD8008"))+
  theme_prism(base_size = 10) +facet_wrap(~C+ cell_type, scale = "free_y") +
  scale_y_continuous(expand = expansion(mult = 0.3)) +
  stat_compare_means(method = "wilcox.test",vjust = -.7, 
                     label.x.npc  = 0.5,  size = 3, fontface = "bold") +
  geom_jitter(alpha = 0.5) 
all_Bs
ggsave("plots/lambda_kappa_ratio_stats.png", all_Bs, w = 5, h = 4.3, dpi = 300, unit = "in")
