source("universal.R")
setRLib()
library(tidyverse)
suppressMessages(library(ggpubr))
library(ggpubr)
library(rstatix)
library(tidyr)
library(ggprism)
library(ggforce)

bcr.heavy <- read_csv("bcr_heavy.csv")
bcr.heavy$V <- gsub("\\-.*", "", bcr.heavy$V)
bcr.heavy$V <- gsub("S.*", "", bcr.heavy$V)
bcr.heavy$V <- gsub("IGHV", "VH", bcr.heavy$V)
bcr.heavy$V <- factor(bcr.heavy$V, levels = c("VH1", "VH2", "VH3", "VH4", "VH5", "VH6", "VH7", "VH8", "VH9", "VH10", "VH11", "VH12", "VH13", "VH14", "VH15", "VH16"))
bcr.heavy$D <- gsub("\\-.*", "", bcr.heavy$D)
bcr.heavy$D <- factor(gsub("IGHD", "DH", bcr.heavy$D))
bcr.heavy$J <- gsub("\\*.*", "", bcr.heavy$J)
bcr.heavy$J <- factor(gsub("IGHJ", "JH", bcr.heavy$J))
bcr.heavy$condition <- factor(bcr.heavy$condition, levels = c("WT", "KO"))
bcr.heavy <- bcr.heavy %>% mutate(cell_type = recode_factor(cell_type, "B1a" = "B-1a","B1b"= "B-1b"))
bcr.h <- bcr.heavy %>% select(count, frequency, CDR3aa, V,D,J, condition, cell_type)
bcr.h <- bcr.h %>% unique
bcr.h_ <- bcr.h %>% group_by(V, J, CDR3aa, condition, cell_type) %>%
  summarize(counts = sum(count)) %>% arrange(-counts)

createDir("./circos_plots")

cell = "B-1a"

B_cell <- bcr.h_ %>% filter(cell_type == cell)
library(circlize)
B_cell <- B_cell %>% ungroup %>% drop_na %>% 
  group_by(condition) %>% mutate(frequency = counts/sum(counts)) %>% ungroup %>% unique()
B_cell$percentages <- B_cell$frequency * 100
###code that might be helpful 
#https://github.com/immunomind/immunarch/issues/103
circos.clear()
lim = 1.2
lims = c(-lim,lim)
f_size = 0.67
.adj = c(-0.8, 0.7)
###ALL VDJ
source("./viz_genus.R")
names(vivid_colors) <- c("JH1","JH2","JH3","JH4","VH12","VH1", "VH2", "VH3", "VH4", "VH5", "VH6", "VH7", 
                         "VH8", "VH9", "VH10", "VH11", "VH12", "VH13", "VH14", "VH15", "VH16")
B_cell_KO <- B_cell %>% filter(condition == "KO")
do_circos(filename = paste0("./circos_plots/", cell, "_KO.pdf"), 
          lims = lims, x = B_cell_KO, font_size = f_size, adj = .adj )
circos.clear()

B_cell_WT <- B_cell %>% filter(condition == "WT")
do_circos(filename =  paste0("./circos_plots/", cell, "_WT.pdf"), 
          lims = lims, x = B_cell_WT, font_size = f_size, adj = .adj )
circos.clear()

#### replicated CDR3

B_cell_KO <- B_cell %>% filter(condition == "KO", frequency > 0.01)
do_circos(filename = paste0("./circos_plots/", cell, "_KO_replicated.pdf"), 
          lims = lims, x = B_cell_KO, font_size = f_size, adj = .adj )
circos.clear()
B_cell_WT <- B_cell %>% filter(condition == "WT", frequency > 0.01)
do_circos(filename = paste0("./circos_plots/", cell, "_WT_replicated.pdf"), 
          lims = lims, x = B_cell_WT, font_size = f_size, adj = .adj )
circos.clear()


#### add quantification? 
library(tidyr)

cell = "B-1b"

  
### Make sure you're working with per thousnd cdr3
B_cell_ <- bcr.h_ %>% filter(cell_type == cell)
B_cell <- B_cell_ %>% ungroup %>% 
  group_by( condition, V, J) %>%
  summarize(counts = sum(counts)) %>% 
  ungroup
B_cell <- B_cell %>% ungroup %>% drop_na %>% 
  group_by(condition) %>% mutate(tot.counts = sum(counts), frequency = counts/sum(counts)) %>% 
  ungroup %>% unique()
B_cell_connections <- B_cell %>%
  tidyr::unite(V_J, V, J, sep = "-", remove = FALSE) %>% 
  filter(J != ".", V != ".") #%>% mutate(V_J = factor(V_J))
connections <- B_cell_ %>% group_by(condition) %>% mutate(tot.counts = sum(counts)) %>% ungroup %>% 
  group_by(CDR3aa, condition) %>% mutate(frequency = counts/tot.counts) %>% 
  filter(frequency > 0.01) %>% tidyr::unite(V_J, V, J, sep = "-", remove = FALSE)  %>% pull(V_J) %>% unique 
tot.con <- B_cell_connections %>% pull(V_J) %>% unique

options(scipen = 999)
B_cell_connections <- B_cell_connections  %>% group_by(condition) %>%
  mutate(scaling = sum(counts)/1000) %>% 
  mutate(norm_counts = counts/scaling) %>% ungroup %>% 
  group_by(condition, V_J) %>%
  mutate(n = round(sum(norm_counts))) %>% unique %>% ungroup %>% group_by(condition, V_J) %>% 
  summarize(normalized_reads= round(sum(norm_counts)), frequency = sum(frequency))  %>% ungroup

tst <- B_cell_connections %>% select(condition, normalized_reads, V_J) %>% 
  pivot_wider(names_from = condition,
                   values_from = normalized_reads, 
                   values_fill = 0) %>% mutate(p.value = 100000)
# results_chi <- tibble()
# for(x in 1:nrow(tst)){
#   dat <- tst[x, ]
#   if(dat$WT == 0 & dat$KO == 0)
#   {
#     dat$p.value <- 1
#     results_chi <- bind_rows(results_chi, dat)
#     next
#   }
#   res <- chi.test(dat$WT, dat$KO)
#   dat$p.value <- res$p.value
#   results_chi <- bind_rows(results_chi, dat)
# }
# lel <- results_chi %>% adjust_pvalue(method =  "fdr", p.col = "p.value")
library(ggpubr)
# chi.test <- function(a, b) {
#   res <- chisq.test(cbind(a, b), correct = T)
#   res$p.value <- res$p.value * length(connections)
#   return(res)
# }
library(ggprism)

# B_cell_connections %>%  filter(V_J %in% connections)  %>% 
#   ggplot(aes(reorder(condition, normalized_reads), normalized_reads, fill = condition)) + 
#   geom_bar(stat = "identity", position = "dodge") + 
#   scale_fill_manual(values = vivid_colors) + 
#   geom_signif(comparisons = list(c("WT","KO")), test = "chi.test", 
#               y_position = max(B_cell_connections$normalized_reads)*1.1, 
# 
#               map_signif_level = c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05, " " = 100000), 
#               textsize = 4, size = 0,
#               tip_length = 0) +
#   facet_wrap(~V_J, nrow = 1) + coord_cartesian(ylim = c(0, max(B_cell_connections$normalized_reads)*1.3)) + 
#   theme_prism(base_size = 11) %+replace% 
#   theme(strip.text.x = element_text(margin = margin(b = 11/3), size = 10)) + 
#   ylab("Normalized counts") + xlab("") 
# ggsave(filename = paste0("~/bcr/circos_plots/", cell, "_quantification.pdf"), dpi=300, w = 10, h = 3, unit = "in")
# ggsave(filename = paste0("~/bcr/circos_plots/", cell, "_quantification.jpeg"), dpi=300, w = 10, h = 3, unit = "in")

names(vivid_colors) <-NULL
B_cell_connections %>%  filter(V_J %in% connections)  %>% 
  ggplot(aes(reorder(condition, normalized_reads), frequency*100, fill = condition)) + 
  geom_bar(stat = "identity", position = "dodge",  color = "black", linewidth = 1) + 
  scale_fill_manual(values = vivid_colors) + 
  theme_prism(base_size = 11) %+replace% 
  theme(strip.text.x = element_text(margin = margin(b = 11/3), size = 10)) + 
  ylab("% abundance") + xlab("") + facet_wrap(~V_J, nrow = 1) 
ggsave(filename = paste0("./circos_plots/", cell, "_quantification.pdf"), dpi=300, w = 10, h = 3, unit = "in")
ggsave(filename = paste0("./circos_plots/", cell, "_quantification.jpeg"), dpi=300, w = 10, h = 3, unit = "in")







