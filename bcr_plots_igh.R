library(tidyverse)
suppressMessages(library(ggpubr))
library(ggpubr)
library(rstatix)
library(tidyr)
library(ggprism)


safe_colorblind_palette <- c("#88CCEE",  "#DDCC77", "#CC6677","#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", 
                             "#E0E678", "#E3B476",  "#FFED61", "#8B8BD9", "#C1F5ED",  "#E6C3F7",
                             "#BF1515")
# files_light <- list.files('old_R/', pattern = '*light*.csv')
# files_heavy <- list.files('old_R/', pattern = '*heavy*.csv')

CompareGroups <- function(dat,col,metric, titile){
  
  p.value <- wilcox.test(get(col) ~ condition, dat)$p.value
  p.value <- sprintf(p.value, fmt = '%#.3f')
  par(oma=c(2,2,2,2))
  y_pos <- max(dat[,col])*1.2
  gp <- ggplot(dat, aes_string(x = factor(dat$condition), y = col, fill = 'condition'))+
    geom_boxplot(lwd=0.3,size=0.3,outlier.size=-1,outlier.shape = NA, width=0.3,alpha=0.7,position = position_dodge2(preserve = "single"))+
    geom_jitter(shape=16, position=position_dodge(0.5), size = 1) +
    ylab(metric)+ 
    xlab("") + 
    scale_fill_manual(values= cls) +
    ggtitle(paste(titile))+
    scale_y_continuous(expand = expansion(mult = 0.2)) +
    theme_prism(base_size =  10) %+replace% theme(plot.title = element_text(size = 10, vjust = 2)) +  guides(fill=FALSE)
  
  if(p.value <0.05)
    gp<- gp+ geom_signif(comparisons = list(c("WT", "KO")), 
                         map_signif_level=TRUE)
  return(gp)
}

file_save <- function(.fname){
  filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/", .fname)
  return(filename)
}
samples <- read_csv("../samples.csv")
rownames(samples) <- samples$library
# filelist.samples.light <- sapply(rownames(samples), function(x) grep(paste0("\\b",x,"_"), files_light, value = TRUE))
# filelist.samples.heavy <- sapply(rownames(samples), function(x) grep(paste0("\\b",x,"_"), files_heavy, value = TRUE))
# 
# merge_process <- function(filelist, sep){
#   output <- list()
#   
#   for( i in names(filelist)){
#     print (filelist[i])
#     file <- filelist[[i]]
#     dat <- read_delim(file, delim = sep)  
#     
#     if(dim(dat)[1] != 0 ) {
#       output[[i]]<- dat  %>%
#         mutate(condition = as.character(samples[i,"condition"]), cell_type = as.character(samples[i,"cell_type"]))
#     } else { output[[i]] <- dat}
#     
#     mat <- do.call(rbind.data.frame, output)
#   }
#   return(mat)
# }
#setwd("old_R/")
# bcr.light <- merge_process(filelist.samples.light, sep = ',')
# bcr.heavy  <- merge_process(filelist.samples.heavy, sep = ',') 
# 
# bcr.heavy$condition <- factor(bcr.heavy$condition, levels = c("WT", "KO"))
# v <- bcr.heavy %>% filter(V != ".") %>% dplyr::group_by(sample, condition, cell_type, V) %>%
#   summarize(freq = sum(frequency)) %>%
#   ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
# tmp1 <- v %>% complete(nesting(smp, condition, cell_type), V) %>% replace(is.na(.), 0)
# 
# cls = c("#0000FF", "#FD8008", "#6db32c", "#cc322a")
# if (!dir.exists("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/")){
#   dir.create("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/")
# }else{
#   print("dir exists")
# }
bcr.heavy <- read_csv("bcr_heavy.csv.gz")

bcr.heavy$condition <- factor(bcr.heavy$condition, levels = c("WT", "KO"))
bcr.heavy <- bcr.heavy %>% mutate(cell_type = recode_factor(cell_type, "B1a" = "B-1a","B1b"= "B-1b"))
#write_csv(bcr.heavy, "../bcr_heavy.csv")


############################# V #############################
#############################################################

v <- bcr.heavy %>% filter(V != ".") %>% dplyr::group_by(sample, condition, cell_type, V) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
tmp1 <- v %>% complete(nesting(smp, condition, cell_type), V) %>% replace(is.na(.), 0)
B1a_p_val <- tmp1 %>% filter(cell_type == "B1a") %>% 
  dplyr::group_by(V) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(V, Wilcox = w$p.value)

sigs <- B1a_p_val %>% filter(Wilcox < 0.05) %>% pull(V)
v_in_B1a <- tmp1 %>% filter(cell_type == "B1a") %>% group_by(condition, V) %>% 
  summarize(avg_freq = mean(freq)) %>% ungroup %>% 
  unique %>% slice(1:15) %>%  pull(V)
B1a_V <- tmp1 %>% filter(cell_type == "B1a", V %in% v_in_B1a) %>% ggplot(aes(x = fct_reorder(V, desc(freq)), y = freq, fill = condition )) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = condition, alpha = 0.5, shape = condition)) +
  scale_fill_manual(values = cls) +
  scale_color_manual(values = cls) +
  theme_prism(base_size = 10) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  ylab("% of all sequences") + xlab("")  +
  ggtitle("V genes") +
  stat_compare_means(label = "p.signif", hide.ns = T, fontface = "bold", size = 6, vjust = 1) + 
  guides(alpha = "none") +
  scale_y_continuous(labels = scales::percent)
B1a_V
ggsave(plot = B1a_V, filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/", "B1a_v_trust4.tiff"),
       height = 3.3, width = 6.6, unit = "in", dpi =300)


B1b_p_val <- tmp1 %>% filter(cell_type == "B1b") %>% 
  dplyr::group_by(V) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(V, Wilcox = w$p.value)


sigs <- B1b_p_val %>% filter(Wilcox < 0.05) %>% pull(V)
v_in_B1b <- tmp1 %>% filter(cell_type == "B1b") %>% group_by(condition, V) %>% 
  summarize(avg_freq = mean(freq)) %>% ungroup %>% 
  arrange(desc(avg_freq)) %>% select(V) %>%
  unique %>% slice(1:15) %>%  pull(V)
B1b_V <- tmp1 %>% filter(cell_type == "B1b", V %in% v_in_B1b) %>% ggplot(aes(x = fct_reorder(V, desc(freq)), y = freq, fill = condition )) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = condition, alpha = 0.5, shape = condition)) +
  scale_fill_manual(values = cls) +
  scale_color_manual(values = cls) +
  theme_prism(base_size = 10) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  ylab("% of all sequences") + xlab("")  +
  ggtitle("V genes") +
  stat_compare_means(label = "p.signif", hide.ns = T, fontface = "bold", size = 6, vjust = 1) + 
  guides(alpha = "none")+
  scale_y_continuous(labels = scales::percent)
B1b_V
ggsave(plot = B1b_V, filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/", "B1b_v_trust4.tiff"),
       height = 3.3, width = 6.6, unit = "in", dpi =300)

############################# D #############################
#############################################################


D <- bcr.heavy %>% filter(V != ".") %>% dplyr::group_by(sample, condition, cell_type, D) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
tmp1 <- D %>% complete(nesting(smp, condition, cell_type), D) %>% replace(is.na(.), 0)

B1a_p_val_d <- tmp1 %>% filter(cell_type == "B1a") %>% 
  dplyr::group_by(D) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(D, Wilcox = w$p.value)


sigs <- B1a_p_val_d %>% filter(Wilcox < 0.05) %>% pull(D)
d_in_B1a <- tmp1 %>% filter(cell_type == "B1a") %>% group_by(condition, D) %>% 
  summarize(avg_freq = mean(freq)) %>% ungroup %>% 
  arrange(desc(avg_freq)) %>% select(D) %>%
  unique %>% slice(1:15) %>%  pull(D)
B1a_D <- tmp1 %>% filter(cell_type == "B1a", D %in% d_in_B1a) %>% ggplot(aes(x = fct_reorder(D, desc(freq)), y = freq, fill = condition )) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = condition, alpha = 0.5, shape = condition)) +
  scale_fill_manual(values = cls) +
  scale_color_manual(values = cls) +
  theme_prism(base_size = 10) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  ylab("% of all sequences") + xlab("")  +
  ggtitle("D genes") +
  stat_compare_means(label = "p.signif", hide.ns = T, fontface = "bold", size = 6, vjust = 1) + 
  guides(alpha = "none")+
  scale_y_continuous(labels = scales::percent)
B1a_D
ggsave(plot = B1a_D, filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/", "B1a_d_trust4.tiff"),
       height = 3.3, width = 6.6, unit = "in", dpi =300)

D <- bcr.heavy %>% filter(V != ".") %>% dplyr::group_by(sample, condition, cell_type, D) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
tmp1 <- D %>% complete(nesting(smp, condition, cell_type), D) %>% replace(is.na(.), 0)

B1b_p_val_d <- tmp1 %>% filter(cell_type == "B1b") %>% 
  dplyr::group_by(D) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(D, Wilcox = w$p.value)
sigs <- B1b_p_val_d %>% filter(Wilcox < 0.05) %>% pull(D)
d_in_B1b <- tmp1 %>% filter(cell_type == "B1b") %>% group_by(condition, D) %>% 
  summarize(avg_freq = mean(freq)) %>% ungroup %>% 
  arrange(desc(avg_freq)) %>% select(D) %>%
  unique %>% slice(1:15) %>%  pull(D)
B1b_D <- tmp1 %>% filter(cell_type == "B1b", D %in% d_in_B1b) %>% ggplot(aes(x = fct_reorder(D, desc(freq)), y = freq, fill = condition )) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = condition, alpha = 0.5, shape = condition)) +
  scale_fill_manual(values = cls) +
  scale_color_manual(values = cls) +
  theme_prism(base_size = 10) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  ylab("% of all sequences") + xlab("")  +
  ggtitle("D genes") +
  stat_compare_means(label = "p.signif", hide.ns = T, fontface = "bold", size = 6, vjust = 1) + 
  guides(alpha = "none")+
  scale_y_continuous(labels = scales::percent)
B1b_D
ggsave(plot = B1b_D, filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/", "B1b_d_trust4.tiff"),
       height = 3.3, width = 6.6, unit = "in", dpi =300)

########### ############## ############## ############## ############## ############## ############## ############## ############## 
############## ############## J ############## ############## ############## ############## ############## ############## ############## 
############## ############## ############## ############## ############## ############## ############## ############## ############## 


J <- bcr.heavy %>% filter(J != ".") %>% dplyr::group_by(sample, condition, cell_type, J) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
tmp1 <- J %>% complete(nesting(smp, condition, cell_type), J) %>% replace(is.na(.), 0)
B1a_p_val_J <- tmp1 %>% filter(cell_type == "B1a") %>% 
  dplyr::group_by(J) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(J, Wilcox = w$p.value)
sigs <- B1a_p_val_J %>% filter(Wilcox < 0.05) %>% pull(J)
J_in_B1a <- tmp1 %>% filter(cell_type == "B1a") %>% group_by(condition, J) %>% 
  summarize(avg_freq = mean(freq)) %>% ungroup %>% 
  arrange(desc(avg_freq)) %>% select(J) %>%
  unique %>% slice(1:15) %>%  pull(J)
B1a_J <- tmp1 %>% filter(cell_type == "B1a", J %in% J_in_B1a) %>% ggplot(aes(x = fct_reorder(J, desc(freq)), y = freq, fill = condition )) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = condition, alpha = 0.5, shape = condition)) +
  scale_fill_manual(values = cls) +
  scale_color_manual(values = cls) +
  theme_prism(base_size = 10) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  ylab("% of all sequences") + xlab("")  +
  ggtitle("J genes") +
  stat_compare_means(label = "p.signif", hide.ns = T, fontface = "bold", size = 6, vjust = 1) + 
  guides(alpha = "none")+
  scale_y_continuous(labels = scales::percent)
B1a_J
ggsave(plot = B1a_J, filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/", "B1a_J_trust4.tiff"),
       height = 3.3, width = 6.6, unit = "in", dpi =300)

B1b_p_val_J <- tmp1 %>% filter(cell_type == "B1b") %>% 
  dplyr::group_by(J) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(J, Wilcox = w$p.value)
sigs <- B1b_p_val_J %>% filter(Wilcox < 0.05) %>% pull(J)
J_in_B1b <- tmp1 %>% filter(cell_type == "B1b") %>% group_by(condition, J) %>% 
  summarize(avg_freq = mean(freq)) %>% ungroup %>% 
  arrange(desc(avg_freq)) %>% select(J) %>%
  unique %>% slice(1:15) %>%  pull(J)
B1b_J <- tmp1 %>% filter(cell_type == "B1b", J %in% J_in_B1b) %>% ggplot(aes(x = fct_reorder(J, desc(freq)), y = freq, fill = condition )) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = condition, alpha = 0.5, shape = condition)) +
  scale_fill_manual(values = cls) +
  scale_color_manual(values = cls) +
  theme_prism(base_size = 10) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  ylab("% of all sequences") + xlab("")  +
  ggtitle("J genes") +
  stat_compare_means(label = "p.signif", hide.ns = T, fontface = "bold", size = 6, vjust = 1) + 
  guides(alpha = "none")+
  scale_y_continuous(labels = scales::percent)
B1b_J
ggsave(plot = B1b_J, filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/", "B1b_J_trust4.tiff"),
       height = 3.3, width = 6.6, unit = "in", dpi =300)

############## ############## ############## ############## ############## ############## ############## ############## ############## 
############## ############## ############## ##############  ISOTYPES AND UNIQUE CDR3 ############## ############## ############## 
############## ############## ############## ############## ############## ############## ############## ############## ############## 
bcr.heavy$Class <- ifelse(str_detect(bcr.heavy$C, "IGHG"), "IGHG", bcr.heavy$C)
st.Ig <- bcr.heavy%>% 
  group_by(sample,condition, cell_type) %>%
  mutate(est_clonal_exp_norm = frequency/sum(frequency)) %>%  
  dplyr::filter(Class != ".", Class != "IGHE", cell_type != "B2") %>%
  group_by(sample,condition,cell_type, Class) %>% 
  dplyr::summarise(Num.Ig = sum(est_clonal_exp_norm))  

barplot_igs <- ggplot(st.Ig,aes(x = condition, y = Num.Ig, fill = Class))+
  geom_bar(stat = "identity",position="fill",alpha = 0.8)+
  theme_bw()+
  labs(y = "Normalized Ig Abundance",fill = "Ig") + facet_wrap(~cell_type) +
  scale_fill_manual(values = c("#6db32c", "#cc322a", "#FD8008", "#0000FF", "#6db32c", "#cc322a")) +
  theme_prism()

ggsave(plot = barplot_igs, filename = file_save("barplot_igs.tiff"), h = 6.5, w = 6.2) 
all_Bs <- st.Ig %>% mutate(Class = factor(Class, levels = c("IGHM", "IGHD"))) %>% dplyr::filter(Class != ".", cell_type != "B2") %>%
  ggplot(aes(x = condition, y = Num.Ig, fill = condition))+ 
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  labs(y = "Normalized Ig Abundance",fill = "Ig")+
  scale_fill_manual(values =  cls)+
  theme_prism(base_size = 10) +facet_wrap(~Class + cell_type, scale = "free_y") +
  scale_y_continuous(expand = expansion(mult = 0.3)) +
  stat_compare_means(method = "wilcox.test", hide.ns = T, vjust = -.7, label.x.npc  = 0.5, label = "p.signif", size = 5, fontface = "bold") + geom_jitter(alpha = 0.5) 

igd <- st.Ig %>% mutate(cell_type = recode_factor(cell_type, "B1a" = "B-1a","B1b"= "B-1b")) %>% 
  dplyr::filter(Class == "IGHD", cell_type != "B2") %>%
  ggplot(aes(x = condition, y = Num.Ig, fill = condition))+ 
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  labs(y = "Normalized Ig Abundance",fill = "Ig")+
  scale_fill_manual(values =  cls)+
  theme_prism(base_size = 10) +facet_wrap(~cell_type) +
  scale_y_continuous(expand = expansion(mult = 0.3)) + ggtitle("IgD") +
  stat_compare_means(method = "wilcox.test", hide.ns = T, vjust = -.7, label.x.npc  = 0.5, label = "p.signif", size = 5, fontface = "bold") + geom_jitter(alpha = 0.5) 
igm <- st.Ig%>% mutate(cell_type = recode_factor(cell_type, "B1a" = "B-1a","B1b"= "B-1b")) %>% dplyr::filter(Class == "IGHM", cell_type != "B2") %>%
  ggplot(aes(x = condition, y = Num.Ig, fill = condition))+ 
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  labs(y = "Normalized Ig Abundance",fill = "Ig")+
  scale_fill_manual(values =  cls)+
  theme_prism(base_size = 10) +facet_wrap(~cell_type) +
  scale_y_continuous(expand = expansion(mult = 0.3)) + ggtitle("IgM") +
  stat_compare_means(method = "wilcox.test", hide.ns = T, vjust = -.7, label.x.npc  = 0.5, label = "p.signif", size = 5, fontface = "bold") + geom_jitter(alpha = 0.5) 
plot <- ggarrange(igm, igd, nrow = 2, common.legend = TRUE, legend = "right")+bgcolor("white") +border("white")
ggsave( plot = plot,
        filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/","all_igs_trust4.tiff"), 
        height = 5, width = 4.12, units = "in", dpi = 300)

all_Bs
ggsave( plot = all_Bs,
        filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/","all_igs_trust4.tiff"), 
        height = 4, width = 4.12, units = "in", dpi = 300)
all_Bs2 <- st.Ig %>% mutate(Class = factor(Class, levels = c("IGHM", "IGHD"))) %>% dplyr::filter(Class != ".", cell_type != "B2") %>%
  ggplot(aes(x = condition, y = Num.Ig, fill = condition))+ 
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  labs(y = "Normalized Ig Abundance",fill = "Ig")+
  scale_fill_manual(values =  cls)+
  theme_prism(base_size = 10) +facet_wrap(~Class + cell_type) +
  scale_y_continuous(expand = expansion(mult = 0.1), labels = scales::percent) +
  stat_compare_means(method = "wilcox.test", hide.ns = T, vjust  = 0.1, label.x.npc  = 0.5, label = "p.signif", size = 5, fontface = "bold") + geom_jitter(alpha = 0.5) 

all_Bs2
ggsave( plot = all_Bs2,
        filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/","all_igs_trust4.3.tiff"), 
        height = 4.6, width = 4.12, units = "in", dpi = 300)

B1a.bcr.heavy <- bcr.heavy %>% filter(cell_type == "B1a")
dat <- aggregate(CDR3aa ~ sample + condition, B1a.bcr.heavy, function(x) length(unique(x))) 
B1a_gc <- CompareGroups(dat,"CDR3aa","Number of Unique CDR3", "B1a")

B1b.bcr.heavy <- bcr.heavy %>% filter(cell_type == "B1b")
dat <- aggregate(CDR3aa ~ sample + condition, B1b.bcr.heavy, function(x) length(unique(x))) 
B1b_gc <- CompareGroups(dat,"CDR3aa","Number of Unique CDR3", "B1b")
B1b_gc <- B1b_gc + ylab("")

B2.bcr.heavy <- bcr.heavy %>% filter(cell_type == "B2")
dat <- aggregate(CDR3aa ~ sample + condition, B2.bcr.heavy, function(x) length(unique(x))) 
B2_gc <- CompareGroups(dat,"CDR3aa","Number of Unique CDR3", "B2")


pl <- ggarrange(B1a_gc, B1b_gc,  nrow = 1, common.legend = TRUE, legend = "bottom")
pl
ggsave(plot = pl, file = file_save("unique_CDR3.jpg"),h = 1.8, w = 3.8, dpi = 300)
.V <- bcr.heavy %>% dplyr::group_by(sample, condition, cell_type, V) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
.D <- bcr.heavy %>% dplyr::group_by(sample, condition, cell_type, D) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
.J <- bcr.heavy %>% dplyr::group_by(sample, condition, cell_type, J) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
write_csv(.V, "V.csv")
write_csv(.D, "D.csv")
write_csv(.J, "J.csv")
write_csv(cdr3, "cdr3.csv")
write_csv(bcr.heavy, "heavy_chain.csv")

############## ############## ############## ############## ############## ############## ############## ############## 
############## CDR3
############## ############## ############## ############## ############## ############## ############## ############## 

cdr3 <- bcr.heavy %>% dplyr::group_by(sample, condition, cell_type, CDR3aa) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
tmp1 <- cdr3 %>% complete(nesting(smp, condition, cell_type), CDR3aa) %>% replace(is.na(.), 0)
B1a_p_val_cdr3 <- tmp1 %>% filter(cell_type == "B1a") %>% 
  dplyr::group_by(CDR3aa) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(CDR3aa, Wilcox = w$p.value)
sigs <- B1a_p_val_cdr3 %>% filter(Wilcox < 0.05, CDR3aa != "out_of_frame") %>% pull(CDR3aa)
cdr3_in_B1a <- tmp1 %>% filter(cell_type == "B1a", CDR3aa != "out_of_frame") %>% group_by(condition, CDR3aa) %>% 
  summarize(avg_freq = mean(freq)) %>% ungroup %>% 
  arrange(desc(avg_freq)) %>% select(CDR3aa) %>%
  unique %>% slice(1:25) %>%  pull(CDR3aa)
tmp1 %>% filter(cell_type == "B1a", CDR3aa != "out_of_frame") %>% 
  arrange(desc(freq)) %>% pull(CDR3aa) %>% unique()

B1a_CDR3aa <- tmp1 %>% 
  filter(cell_type == "B1a", CDR3aa %in% cdr3_in_B1a,  CDR3aa != "out_of_frame") %>% 
  arrange(desc(freq)) %>%
  ggplot(aes(x = reorder(CDR3aa, freq), y = freq, fill = condition )) +
  geom_boxplot(outlier.shape = NA) +
  #geom_jitter(aes(color = condition, alpha = 0.5, shape = condition)) +
  scale_fill_manual(values = cls) +
  scale_color_manual(values = cls) +
  theme_prism(base_size = 10) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  ylab("% of all sequences") + xlab("")  +
  coord_flip(ylim = c(0, 0.4)) +
  stat_compare_means(label = "p.signif", hide.ns = T, fontface = "bold", size = 6, label.y = .4, vjust = 0.7) + 
  guides(alpha = "none")+
  scale_y_continuous(labels = scales::percent)

B1a_CDR3aa

ggsave(plot = B1a_CDR3aa, 
       filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/",
                         "B1a_CDR3aa_trust4.tiff"),
       height = 5.84, width = 6.84, unit = "in", dpi =300)

cdr3 <- bcr.heavy %>% dplyr::group_by(sample, condition, cell_type, CDR3aa) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
tmp1 <- cdr3 %>% complete(nesting(smp, condition, cell_type), CDR3aa) %>% replace(is.na(.), 0)
B1b_p_val_cdr3 <- tmp1 %>% filter(cell_type == "B1b") %>% 
  dplyr::group_by(CDR3aa) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(CDR3aa, Wilcox = w$p.value)
sigs <- B1b_p_val_cdr3 %>% filter(Wilcox < 0.05, CDR3aa != "out_of_frame") %>% pull(CDR3aa)
cdr3_in_B1b <- tmp1 %>% filter(cell_type == "B1b", CDR3aa != "out_of_frame") %>% group_by(condition, CDR3aa) %>% 
  summarize(avg_freq = mean(freq)) %>% ungroup %>% 
  arrange(desc(avg_freq)) %>% select(CDR3aa) %>%
  unique %>% slice(1:25) %>%  pull(CDR3aa)
tmp1 %>% filter(cell_type == "B1b", CDR3aa != "out_of_frame") %>% 
  arrange(desc(freq)) %>% pull(CDR3aa) %>% unique()

B1b_CDR3aa <- tmp1 %>% 
  filter(cell_type == "B1b", CDR3aa %in% cdr3_in_B1b,  CDR3aa != "out_of_frame") %>% 
  arrange(desc(freq)) %>%
  ggplot(aes(x = reorder(CDR3aa, freq), y = freq, fill = condition )) +
  geom_boxplot(outlier.shape = NA) +
  #geom_jitter(aes(color = condition, alpha = 0.5, shape = condition)) +
  scale_fill_manual(values = cls) +
  scale_color_manual(values = cls) +
  theme_prism(base_size = 10) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  ylab("% of all sequences") + xlab("")  +
  coord_flip(ylim = c(0, 0.4)) +
  stat_compare_means(label = "p.signif", hide.ns = T, fontface = "bold", size = 6, label.y = .4, vjust = 0.7) + 
  guides(alpha = "none")+
  scale_y_continuous(labels = scales::percent)
#+
  #scale_y_continuous(labels = scales::percent)

B1b_CDR3aa
ggsave(plot = B1b_CDR3aa, 
       filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/",
                         "B1b_CDR3aa_trust4.tiff"),
       height = 5.84, width = 6.84, unit = "in", dpi =300)
###############################
pl1 <- bcr.heavy %>% 
  filter(cell_type == "B1a",CDR3aa != "out_of_frame", condition == "WT")  %>%  slice(1:20) %>% 
  ggplot(aes(reorder(CDR3aa, frequency), frequency, fill = frequency)) + 
  geom_bar(stat = "identity") +
  coord_flip() + 
  scale_fill_gradient(high="orange", low="blue") +  coord_flip() +
  scale_y_continuous(labels = scales::percent, limits = c(0, .4)) + 
  xlab("CDR3aa")  + 
  ylab("% of all CDR3") + theme_prism(base_size = 8)

pl2 <- bcr.heavy %>% 
  filter(cell_type == "B1a",CDR3aa != "out_of_frame", condition == "KO")  %>%  slice(1:20) %>% 
  ggplot(aes(reorder(CDR3aa, frequency), frequency, fill = frequency)) + 
  geom_bar(stat = "identity") +
  coord_flip() + 
  scale_fill_gradient(high="orange", low="blue") +  coord_flip() +
  scale_y_continuous(labels = scales::percent, limits = c(0, .4)) + 
  xlab("CDR3aa")  + 
  ylab("% of all CDR3") + theme_prism(base_size = 8)

plot <- ggarrange(pl1, pl2, nrow = 1, common.legend = TRUE, legend = "right")+bgcolor("white") +border("white")
ggsave(file_save("CDR3aa_B1a.png"), plot = plot,  height = 4.5, width = 8, unit = "in")

pl1 <- bcr.heavy %>% 
  filter(cell_type == "B1b",CDR3aa != "out_of_frame", condition == "WT")  %>%  slice(1:20) %>% 
  ggplot(aes(reorder(CDR3aa, frequency), frequency, fill = frequency)) + 
  geom_bar(stat = "identity") +
  coord_flip() + 
  scale_fill_gradient(high="orange", low="blue") +  coord_flip() +
  scale_y_continuous(labels = scales::percent, limits = c(0, .4)) + 
  xlab("CDR3aa")  + 
  ylab("% of all CDR3") + theme_prism(base_size = 8)

pl2 <- bcr.heavy %>% 
  filter(cell_type == "B1b",CDR3aa != "out_of_frame", condition == "KO")  %>%  slice(1:20) %>% 
  ggplot(aes(reorder(CDR3aa, frequency), frequency, fill = frequency)) + 
  geom_bar(stat = "identity") +
  coord_flip() + 
  scale_fill_gradient(high="orange", low="blue") +  coord_flip() +
  scale_y_continuous(labels = scales::percent, limits = c(0, .4)) + 
  xlab("CDR3aa")  + 
  ylab("% of all CDR3") + theme_prism(base_size = 8)

plot <- ggarrange(pl1, pl2, nrow = 1, common.legend = TRUE, legend = "right")+bgcolor("white") +border("white")
ggsave(file_save("CDR3aa_B1b.png"), plot = plot,  height = 4.5, width = 8, unit = "in")

############## ############## ############## ############## ############## ############## ############## ############## ############## 
############## ############## ##############  IGM
############## ############## ############## ############## ############## ############## ############## ############## ############## 
igms <- bcr.heavy %>% filter(cell_type == "B1a", Class == "IGHM")
igms <- igms %>% dplyr::group_by(sample, condition, cell_type, CDR3aa) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
tmp1 <- igms %>% complete(nesting(smp, condition, cell_type), CDR3aa) %>% replace(is.na(.), 0)
igms_p_val_cdr3 <- tmp1 %>% filter(cell_type == "B1a") %>% 
  dplyr::group_by(CDR3aa) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(CDR3aa, Wilcox = w$p.value)
sigs <- igms_p_val_cdr3 %>% filter(Wilcox < 0.05, CDR3aa != "out_of_frame") %>% pull(CDR3aa)
igms_in_B1b <- tmp1 %>% filter(cell_type == "B1a", CDR3aa != "out_of_frame") %>% group_by(condition, CDR3aa) %>% 
  summarize(avg_freq = mean(freq)) %>% ungroup %>% 
  arrange(desc(avg_freq)) %>% select(CDR3aa) %>%
  unique %>% slice(1:25) %>%  pull(CDR3aa)
tmp1 %>% filter(cell_type == "B1a", CDR3aa != "out_of_frame") %>% 
  arrange(desc(freq)) %>% pull(CDR3aa) %>% unique()


igms <- tmp1 %>% 
  filter(cell_type == "B1a", CDR3aa %in% igms_in_B1b,  CDR3aa != "out_of_frame") %>% 
  arrange(desc(freq)) %>%
  ggplot(aes(x = reorder(CDR3aa, freq), y = freq, fill = condition )) +
  geom_boxplot(outlier.shape = NA) +
  #geom_jitter(aes(color = condition, alpha = 0.5, shape = condition)) +
  scale_fill_manual(values = cls) +
  scale_color_manual(values = cls) +
  theme_prism(base_size = 10) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  ylab("% of all sequences") + xlab("")  +
  coord_flip(ylim = c(0, 1)) +
  stat_compare_means(label = "p.signif", hide.ns = T, fontface = "bold", size = 6, label.y = .4, vjust = 0.7) + 
  guides(alpha = "none")+
  scale_y_continuous(labels = scales::percent)
igms


###### ALternative V analysis
all_V <- bcr.heavy
all_V$V <- gsub("\\-.*", "", bcr.heavy$V)
all_V$V <- gsub("S.*", "", all_V$V)
all_V$V <- gsub("IGHV", "VH", all_V$V)
v <- all_V %>% filter(V != ".") %>% dplyr::group_by(sample, condition, cell_type, V) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
v$V <- factor(v$V, levels = c("VH1", "VH2", "VH3", "VH4", "VH5", "VH6", "VH7", "VH8", "VH9", "VH10", "VH11", "VH12", "VH13", "VH14", "VH15", "VH16"))
tmp1 <- v %>% complete(nesting(smp, condition, cell_type), V) %>% replace(is.na(.), 0)
B1a_p_val <- tmp1 %>% filter(cell_type == "B1a") %>% 
  dplyr::group_by(V) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(V, Wilcox = w$p.value)

sigs <- B1a_p_val %>% filter(Wilcox < 0.05) %>% pull(V)
B1a_V <- tmp1 %>% filter(cell_type == "B1a") %>% ggplot(aes(x = V, y = freq, fill = condition )) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = condition, alpha = 0.5, shape = condition)) +
  scale_fill_manual(values = cls) +
  scale_color_manual(values = cls) +
  theme_prism(base_size = 10) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  ylab("% of all sequences") + xlab("")  +
  ggtitle("V genes") +
  stat_compare_means(label = "p.signif", hide.ns = T, fontface = "bold", size = 6, vjust = 1) + 
  guides(alpha = "none") +
  scale_y_continuous(labels = scales::percent)
B1a_V
ggsave(plot = B1a_V, filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/", "B1a_v_VH.tiff"),
       height = 3.3, width = 6.6, unit = "in", dpi =300)


all_V <- bcr.heavy
all_V$V <- gsub("\\-.*", "", bcr.heavy$V)
all_V$V <- gsub("S.*", "", all_V$V)
all_V$V <- gsub("IGHV", "VH", all_V$V)
v <- all_V %>% filter(V != ".") %>% dplyr::group_by(sample, condition, cell_type, V) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
v$V <- factor(v$V, levels = c("VH1", "VH2", "VH3", "VH4", "VH5", "VH6", "VH7", "VH8", "VH9", "VH10", "VH11", "VH12", "VH13", "VH14", "VH15", "VH16"))
tmp1 <- v %>% complete(nesting(smp, condition, cell_type), V) %>% replace(is.na(.), 0)
B1b_p_val <- tmp1 %>% filter(cell_type == "B1b") %>% 
  dplyr::group_by(V) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(V, Wilcox = w$p.value)

sigs <- B1b_p_val %>% filter(Wilcox < 0.05) %>% pull(V)
B1b_V <- tmp1 %>% filter(cell_type == "B1b") %>% ggplot(aes(x = V, y = freq, fill = condition )) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = condition, alpha = 0.5, shape = condition)) +
  scale_fill_manual(values = cls) +
  scale_color_manual(values = cls) +
  theme_prism(base_size = 10) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  ylab("% of all sequences") + xlab("")  +
  ggtitle("V genes") +
  stat_compare_means(label = "p.signif", hide.ns = T, fontface = "bold", size = 6, vjust = 1) + 
  guides(alpha = "none") +
  scale_y_continuous(labels = scales::percent)
B1b_V
ggsave(plot = B1b_V, filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/", "B1b_v_VH.tiff"),
       height = 3.3, width = 6.6, unit = "in", dpi =300)


##D
all_D <- bcr.heavy
all_D$D <- gsub("\\-.*", "", bcr.heavy$D)
all_D$D <- gsub("IGHD", "DH", all_D$D)
D <- all_D %>% filter(D != ".") %>% dplyr::group_by(sample, condition, cell_type, D) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
tmp1 <- D %>% complete(nesting(smp, condition, cell_type), D) %>% replace(is.na(.), 0)
B1a_p_val <- tmp1 %>% filter(cell_type == "B1a") %>% 
  dplyr::group_by(D) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(D, Wilcox = w$p.value)

sigs <- B1a_p_val %>% filter(Wilcox < 0.05) %>% pull(D)
B1a_D <- tmp1 %>% filter(cell_type == "B1a") %>% ggplot(aes(x = D, y = freq, fill = condition )) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = condition, alpha = 0.5, shape = condition)) +
  scale_fill_manual(values = cls) +
  scale_color_manual(values = cls) +
  theme_prism(base_size = 10) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  ylab("% of all sequences") + xlab("")  +
  ggtitle("D genes") +
  stat_compare_means(label = "p.signif", hide.ns = T, fontface = "bold", size = 6, vjust = 1) + 
  guides(alpha = "none") +
  scale_y_continuous(labels = scales::percent)
B1a_D
ggsave(plot = B1a_D, filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/", "B1a_D_DH.tiff"),
       height = 3.3, width = 6.6, unit = "in", dpi =300)


B1b_p_val <- tmp1 %>% filter(cell_type == "B1b") %>% 
  dplyr::group_by(D) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(D, Wilcox = w$p.value)

sigs <- B1b_p_val %>% filter(Wilcox < 0.05) %>% pull(D)
B1b_D <- tmp1 %>% filter(cell_type == "B1b") %>% ggplot(aes(x = D, y = freq, fill = condition )) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = condition, alpha = 0.5, shape = condition)) +
  scale_fill_manual(values = cls) +
  scale_color_manual(values = cls) +
  theme_prism(base_size = 10) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  ylab("% of all sequences") + xlab("")  +
  ggtitle("D genes") +
  stat_compare_means(label = "p.signif", hide.ns = T, fontface = "bold", size = 6, vjust = 1) + 
  guides(alpha = "none") +
  scale_y_continuous(labels = scales::percent)
B1b_D
ggsave(plot = B1b_D, filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/", "B1b_D_DH.tiff"),
       height = 3.3, width = 6.6, unit = "in", dpi =300)

#### J
all_J <- bcr.heavy
all_J$J <- gsub("\\*.*", "", bcr.heavy$J)
all_J$J <- gsub("IGHJ", "JH", all_J$J)
J <- all_J %>% filter(J != ".") %>% dplyr::group_by(sample, condition, cell_type, J) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
tmp1 <- J %>% complete(nesting(smp, condition, cell_type), J) %>% replace(is.na(.), 0)
B1a_p_val <- tmp1 %>% filter(cell_type == "B1a") %>% 
  dplyr::group_by(J) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(J, Wilcox = w$p.value)

sigs <- B1a_p_val %>% filter(Wilcox < 0.05) %>% pull(J)
B1a_J <- tmp1 %>% filter(cell_type == "B1a") %>% ggplot(aes(x = J, y = freq, fill = condition )) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = condition, alpha = 0.5, shape = condition)) +
  scale_fill_manual(values = cls) +
  scale_color_manual(values = cls) +
  theme_prism(base_size = 10) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  ylab("% of all sequences") + xlab("")  +
  ggtitle("J genes") +
  stat_compare_means(label = "p.signif", hide.ns = T, fontface = "bold", size = 6, vjust = 1) + 
  guides(alpha = "none") +
  scale_y_continuous(labels = scales::percent)
B1a_J
ggsave(plot = B1a_J, filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/", "B1a_J_JH.tiff"),
       height = 3.3, width = 6.6, unit = "in", dpi =300)


B1b_p_val <- tmp1 %>% filter(cell_type == "B1b") %>% 
  dplyr::group_by(J) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(J, Wilcox = w$p.value)

sigs <- B1b_p_val %>% filter(Wilcox < 0.05) %>% pull(J)
B1b_J <- tmp1 %>% filter(cell_type == "B1b") %>% ggplot(aes(x = J, y = freq, fill = condition )) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color = condition, alpha = 0.5, shape = condition)) +
  scale_fill_manual(values = cls) +
  scale_color_manual(values = cls) +
  theme_prism(base_size = 10) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  ylab("% of all sequences") + xlab("")  +
  ggtitle("J genes") +
  stat_compare_means(label = "p.signif", hide.ns = T, fontface = "bold", size = 6, vjust = 1) + 
  guides(alpha = "none") +
  scale_y_continuous(labels = scales::percent)
B1b_J
ggsave(plot = B1b_J, filename = paste0("/project/mcnamara-lab/Maria/dennis_rnaseq/bcr_take2_trust4/figures/", "B1b_J_JH.tiff"),
       height = 3.3, width = 6.6, unit = "in", dpi =300)




###########
################################ V presence in replicated sequences
##########################
colo <- c('#4363d8', '#f58231','#e6194b', '#3cb44b', '#ffe119',  '#911eb4',
          '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
          '#9a6324', '#fffac8', '#800000', '#aaffc3',
          '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
bcr.heavy.v <- bcr.heavy
bcr.heavy.v$V <- gsub("\\-.*", "", bcr.heavy.v$V)
bcr.heavy.v$V <- gsub("S.*", "", bcr.heavy.v$V)
bcr.heavy.v$V <- gsub("IGHV", "VH", bcr.heavy.v$V)

CDR3aa <- bcr.heavy.v %>% filter(cell_type == "B-1a") %>%
  select(cell_type, condition,CDR3aa, CDR3nt, count) %>% 
  group_by(cell_type, condition) %>%
  mutate(seq_n = sum(count)) %>%
  ungroup %>% 
  group_by(cell_type, condition, CDR3aa)  %>% 
  summarize(c_CDR3 = sum(count), f_CDR3 = sum(count)/seq_n ) %>%
  ungroup %>% distinct %>% arrange(-f_CDR3)

CDR3aa.V <- bcr.heavy.v %>% filter(cell_type == "B-1a") %>%
  select(cell_type, condition,CDR3aa, CDR3nt, count, V) %>% 
  group_by(cell_type, condition) %>%
  mutate(seq_n = sum(count)) %>%
  ungroup %>% 
  group_by(cell_type, condition, CDR3aa, V)  %>% 
  summarize(c_CDR3 = sum(count), f_CDR3 = sum(count)/seq_n ) %>%
  ungroup %>% distinct %>% arrange(-f_CDR3) %>% ungroup
  

CDR3_props <- CDR3aa.V %>% select(CDR3aa, condition, f_CDR3, c_CDR3) %>% group_by(condition) %>%  unique %>% ungroup
table(CDR3_props$condition, CDR3_props$f_CDR3 > 0.01)
seq_clones_WT <- sum((CDR3_props[CDR3_props$f_CDR3 > 0.01 & CDR3_props$condition == "WT",])$f_CDR3)
seq_clones_KO <- sum((CDR3_props[CDR3_props$f_CDR3 > 0.01 & CDR3_props$condition == "KO",])$f_CDR3)
CDR3_props$CDR3 <- ifelse(CDR3_props$f_CDR3 > 0.01, CDR3_props$CDR3aa, "not replicated")
CDR3_props <- CDR3_props %>%  group_by(condition, CDR3) %>%  summarize(f = sum(f_CDR3), c = sum(c_CDR3))



cd1 <- CDR3_props %>% filter(condition == "WT") %>% ggplot(aes(x = "", y=c,fill = reorder(CDR3, -c))) +
  geom_bar(stat="identity", width=0.001, color="white") + labs(fill = "CDR3 AA sequence (B-1a WT)")+
  coord_polar("y", start=0) + 
  scale_fill_manual(values = colo) +
  theme_classic(base_size = 14) %+replace% theme(line = element_blank(), axis.text = element_blank()) +
  ylab("") + xlab("") 
ggsave(plot = cd1, filename = file_save("piechart_B1a_WT_cdr3.tiff"), dpi = 300, w = 7, h = 4.5)
cd2 <- CDR3_props %>% filter(condition == "KO") %>% ggplot(aes(x = "", y=c, fill = reorder(CDR3, -c))) +
  geom_bar(stat="identity", width=0.001, color="white") + labs(fill = "CDR3 AA sequence (B-1a KO)")+
  coord_polar("y", start=0) +
  scale_fill_manual(values = colo) +
  theme_classic(base_size = 14) %+replace% theme(line = element_blank(), axis.text = element_blank()) +
  ylab("") + xlab("") 
ggsave(plot = cd2, filename = file_save("piechart_B1a_KO_cdr3.tiff"), dpi = 300, w = 7, h = 4.5)

CDR3_props <- CDR3aa.V %>% select(CDR3aa, condition, f_CDR3, c_CDR3) %>% group_by(condition) %>%  unique %>% ungroup
CDR3_props$rep<- ifelse(CDR3_props$f_CDR3 > 0.01, "replicated", "not replicated")
f_test <- CDR3_props %>% mutate(f = round(f_CDR3 * 100))
cond1 <- rep(f_test$condition ,f_test$f)
cond2 <- rep(f_test$rep ,f_test$f)
cond <- tibble(condition = cond1, rep = cond2)
table(cond$condition, cond$rep)
chisq.test(cond$condition, cond$rep)
table(CDR3_props$condition, CDR3_props$rep, CDR3_props$c_CDR3)
c_test <- CDR3_props %>% mutate(c = round(c_CDR3))

cond1 <- rep(c_test$condition ,c_test$c)
cond2 <- rep(c_test$rep ,c_test$c)
cond <- tibble(condition = cond1, rep = cond2)
table(cond$condition, cond$rep)
rep_test <-CDR3_props  %>% group_by(condition) %>% mutate(total = n()) %>% ungroup %>% 
  group_by(condition, rep, total) %>%  summarize(n = n()) %>% mutate(perc = n/total)

for(sf in seq(1, 10000000, 100)){
  rep_test_all <-CDR3_props  %>% group_by(condition) %>% mutate(scaling = sum(c_CDR3)/sf) %>% 
    mutate(norm_counts = c_CDR3/scaling) %>% ungroup %>% 
    group_by(condition, rep) %>%
    summarize(n = round(sum(norm_counts))) %>% unique
  test <- matrix(c(rep_test_all$n[1],rep_test_all$n[3],rep_test_all$n[2], rep_test_all$n[4]), nrow = 2)
  resu <- chisq.test(test)
  res <- tibble(scaling_factor = sf, "-log10(p.value)" = -log10(resu$p.value))
  scaling_factors <- bind_rows(scaling_factors, res)
}

sfs <- c(10, 50, 100, 500,
  1000, 5000, 10000,
  50000, 100000, 500000,
  1000000, 5000000, 10000000,
  50000000, 100000000)
library(foreach)
n.cores <- parallel::detectCores() - 1
library(parallel)
options(scipen = 100000000000000000)
cl <- makeCluster(n.cores)

doParallel::registerDoParallel(cl = cl)
foreach::getDoParRegistered()


foreach(sf = seq(1, 100000000, 5000)) %dopar% {
  rep_test_all <-CDR3_props  %>% group_by(condition) %>% mutate(scaling = sum(c_CDR3)/sf) %>% 
    mutate(norm_counts = c_CDR3/scaling) %>% ungroup %>% 
    group_by(condition, rep) %>%
    summarize(n = round(sum(norm_counts))) %>% unique
  test <- matrix(c(rep_test_all$n[1],rep_test_all$n[3],rep_test_all$n[2], rep_test_all$n[4]), nrow = 2)
  resu <- chisq.test(test)
  res <- tibble(scaling_factor = sf, "-log10(p.value)" = -log10(resu$p.value))
  scaling_factors <- bind_rows(scaling_factors, res)
}

chisq.test( CDR3_props$condition, CDR3_props$rep, p = CDR3_props$f_CDR3, rescale.p = T)
chisq.test( CDR3_props$condition, CDR3_props$rep)
fisher.test( table(CDR3_props$condition, CDR3_props$rep))
fisher.test( table(CDR3_props$condition, CDR3_props$rep))

CDR3aa <- bcr.heavy.v %>% filter(cell_type == "B-1b") %>%
  select(cell_type, condition,CDR3aa, CDR3nt, count) %>% 
  group_by(cell_type, condition) %>%
  mutate(seq_n = sum(count)) %>%
  ungroup %>% 
  group_by(cell_type, condition, CDR3aa)  %>% 
  summarize(c_CDR3 = sum(count), f_CDR3 = sum(count)/seq_n ) %>%
  ungroup %>% distinct %>% arrange(-f_CDR3)

CDR3aa.V <- bcr.heavy.v %>% filter(cell_type == "B-1b") %>%
  select(cell_type, condition,CDR3aa, CDR3nt, count, V) %>% 
  group_by(cell_type, condition) %>%
  mutate(seq_n = sum(count)) %>%
  ungroup %>% 
  group_by(cell_type, condition, CDR3aa, V)  %>% 
  summarize(c_CDR3 = sum(count), f_CDR3 = sum(count)/seq_n ) %>%
  ungroup %>% distinct %>% arrange(-f_CDR3) %>% ungroup


CDR3_props <- CDR3aa.V %>% select(CDR3aa, condition, f_CDR3, c_CDR3) %>% group_by(condition) %>%  unique %>% ungroup
table(CDR3_props$condition, CDR3_props$f_CDR3 > 0.01)
seq_clones_WT <- sum((CDR3_props[CDR3_props$f_CDR3 > 0.01 & CDR3_props$condition == "WT",])$f_CDR3)
seq_clones_KO <- sum((CDR3_props[CDR3_props$f_CDR3 > 0.01 & CDR3_props$condition == "KO",])$f_CDR3)
CDR3_props$CDR3 <- ifelse(CDR3_props$f_CDR3 > 0.01, CDR3_props$CDR3aa, "not replicated")
CDR3_props <- CDR3_props %>%  group_by(condition, CDR3) %>%  summarize(f = sum(f_CDR3), c = sum(c_CDR3))




cd1 <- CDR3_props %>% filter(condition == "WT") %>% ggplot(aes(x = "", y=c,fill = reorder(CDR3, -c))) +
  geom_bar(stat="identity", width=0.001, color="white") + labs(fill = "CDR3 AA sequence (B-1b WT)")+
  coord_polar("y", start=0) + 
  scale_fill_manual(values = colo) +
  theme_classic(base_size = 14) %+replace% theme(line = element_blank(), axis.text = element_blank()) +
  ylab("") + xlab("") 
ggsave(plot = cd1, filename = file_save("piechart_B1b_WT_cdr3.tiff"), dpi = 300, w = 7, h = 4.5)
cd2 <- CDR3_props %>% filter(condition == "KO") %>% ggplot(aes(x = "", y=c, fill = reorder(CDR3, -c))) +
  geom_bar(stat="identity", width=0.001, color="white") + labs(fill = "CDR3 AA sequence (B-1b KO)")+
  coord_polar("y", start=0) +
  scale_fill_manual(values = colo) +
  theme_classic(base_size = 14) %+replace% theme(line = element_blank(), axis.text = element_blank()) +
  ylab("") + xlab("") 
ggsave(plot = cd2, filename = file_save("piechart_B1b_KO_cdr3.tiff"), dpi = 300, w = 7, h = 4.5)


CDR3_props <- CDR3aa.V %>% select(CDR3aa, condition, f_CDR3, c_CDR3) %>% group_by(condition) %>%  unique %>% ungroup
CDR3_props$rep<- ifelse(CDR3_props$f_CDR3 > 0.01, "replicated", "not replicated")
f_test <- CDR3_props %>% mutate(f = round(f_CDR3 * 100))
cond1 <- rep(f_test$condition ,f_test$f)
cond2 <- rep(f_test$rep ,f_test$f)
cond <- tibble(condition = cond1, rep = cond2)
table(cond$condition, cond$rep)
chisq.test(cond$condition, cond$rep)
table(CDR3_props$condition, CDR3_props$rep, CDR3_props$c_CDR3)
c_test <- CDR3_props %>% mutate(c = round(c_CDR3))

cond1 <- rep(c_test$condition ,c_test$c)
cond2 <- rep(c_test$rep ,c_test$c)
cond <- tibble(condition = cond1, rep = cond2)
table(cond$condition, cond$rep)
rep_test <-CDR3_props  %>% group_by(condition) %>% mutate(total = n()) %>% ungroup %>% 
  group_by(condition, rep, total) %>%  summarize(n = n()) %>% mutate(perc = n/total)
sf = 1000000
rep_test_all <-CDR3_props  %>% group_by(condition) %>% mutate(scaling = sum(c_CDR3)/sf) %>% 
  mutate(norm_counts = c_CDR3/scaling) %>% ungroup %>% 
  group_by(condition, rep) %>%
  summarize(n = round(sum(norm_counts))) %>% unique


fisher.test( table(CDR3_props$condition, CDR3_props$rep))

chisq.test( CDR3_props$condition, CDR3_props$rep, p = CDR3_props$f_CDR3)



####test 
scaling_factors <- tibble()

sf_determin <- function(sf = 1){
  if(sf == 0) return(NA)
  rep_test_all <-CDR3_props  %>% group_by(condition) %>% mutate(scaling = sum(c_CDR3)/sf) %>% 
    mutate(norm_counts = c_CDR3/scaling) %>% ungroup %>% 
    group_by(condition, rep) %>%
    summarize(n = round(sum(norm_counts))) %>% unique
  test <- matrix(c(rep_test_all$n[1],rep_test_all$n[3],rep_test_all$n[2], rep_test_all$n[4]), nrow = 2)
  resu <- chisq.test(test)
  res <- c(scaling_factor = sf, "-log10(p.value)" = -log10(resu$p.value))
  }
res <- lapply(seq(0, 100000, 500), sf_determin)
res.copy <- res
df_res <- as.data.frame(do.call(rbind, res))
x <- sf_determin(2)
df_res <- df_res %>% rbind(x)

x <- sf_determin(10)
df_res <- df_res %>% rbind(x)

x <- sf_determin(100)
df_res <- df_res %>% rbind(x)

x <- sf_determin(200)
df_res <- df_res %>% rbind(x)

df_res <- df_res %>% as_tibble %>% drop_na() %>% mutate(pvalue = `-log10(p.value)`)
df_res %>% ggplot(aes(x = scaling_factor, y = pvalue)) + geom_line() + theme_bw() + ylab("-log10(pvalue)")


for(sf in seq(1, 10000000, 500)){
  rep_test_all <-CDR3_props  %>% group_by(condition) %>% mutate(scaling = sum(c_CDR3)/sf) %>% 
    mutate(norm_counts = c_CDR3/scaling) %>% ungroup %>% 
    group_by(condition, rep) %>%
    summarize(n = round(sum(norm_counts))) %>% unique
  test <- matrix(c(rep_test_all$n[1],rep_test_all$n[3],rep_test_all$n[2], rep_test_all$n[4]), nrow = 2)
  resu <- chisq.test(test)
  res <- tibble(scaling_factor = sf, "-log10(p.value)" = -log10(resu$p.value))
  scaling_factors <- bind_rows(scaling_factors, res)
}

