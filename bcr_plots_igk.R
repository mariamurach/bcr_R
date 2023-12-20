setwd("~/bcr/")
rver <- getRversion()
libdir <- paste0("/home/mm5jy/R/", rver)
if (!dir.exists(libdir)){
  dir.create(libdir)
}else{
  print("dir exists")
}
.libPaths(libdir)

library(tidyverse)
suppressMessages(library(ggpubr))
library(ggpubr)
library(rstatix)
library(tidyr)
library(ggprism)
theme_mm <- function() {
  font <- "Georgia"   #assign font family up front
  
  theme_bw() %+replace%    #replace elements we want to change
    theme(
      panel.background = element_rect(fill = "white",
                                      colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      text = element_text(size = 8),
      plot.title = element_text(size = 8, margin = margin(rep(2,4))),
      axis.line = element_blank(),
      legend.key = element_blank(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 0.5
      ),
      strip.background = element_rect(fill = "white"),
      complete = TRUE,
      axis.title.x = element_blank(),
      
      legend.position="right"
    )
}

theme_mm_lines  <- function() {
  font <- "Georgia"   #assign font family up front
  
  theme_bw() %+replace%    #replace elements we want to change
    theme(
      panel.background = element_rect(fill = "white",
                                      colour = NA),
      panel.grid = element_line(colour = "grey95"),
      text = element_text(size = 8),
      plot.title = element_text(size = 8, margin = margin(rep(2,4))),
      axis.line = element_blank(),
      legend.key = element_blank(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 0.5
      ),
      complete = TRUE,
      axis.title.x = element_blank(),
      strip.background = element_rect(fill = "white",
                                      colour = "grey20"), 
      legend.position="right"
    )
}

safe_colorblind_palette <- c("#88CCEE",  "#DDCC77", "#CC6677","#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", 
                             "#E0E678", "#E3B476",  "#FFED61", "#8B8BD9", "#C1F5ED",  "#E6C3F7",
                             "#BF1515")
files_light <- list.files('old_R/', pattern = '*light*.csv')

CompareGroups <- function(dat,col,metric, titile){
  
  p.value <- wilcox.test(get(col) ~ condition, dat)$p.value
  p.value <- sprintf(p.value, fmt = '%#.3f')
  print(p.value)
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
  filename = paste0("~/bcr/figures_light/", .fname)
  return(filename)
}

if (!dir.exists("~/bcr/figures_light/")){
  dir.create("~/bcr/figures_light/")
}else{
  print("dir exists")
}
samples <- read_csv("samples.csv")
rownames(samples) <- samples$library
filelist.samples.light <- sapply(rownames(samples), function(x) grep(paste0("\\b",x,"_"), files_light, value = TRUE))

merge_process <- function(filelist, sep){
  output <- list()

  for( i in names(filelist)){
    print (filelist[i])
    file <- filelist[[i]]
    dat <- read_delim(file, delim = sep)

    if(dim(dat)[1] != 0 ) {
      output[[i]]<- dat  %>%
        mutate(condition = as.character(samples[i,"condition"]), cell_type = as.character(samples[i,"cell_type"]))
    } else { output[[i]] <- dat}

    mat <- do.call(rbind.data.frame, output)
  }
  return(mat)
}
setwd("~/bcr/old_R/")
bcr.light <- merge_process(filelist.samples.light, sep = ',')
 
cls = c("#0000FF", "#FD8008", "#6db32c", "#cc322a")
if (!dir.exists("~/bcr/figures/")){
  dir.create("~/bcr/figures/")
}else{
  print("dir exists")
}

bcr.light$condition <- factor(bcr.light$condition, levels = c("WT", "KO"))
bcr.light <- bcr.light %>% mutate(cell_type = recode_factor(cell_type, "B1a" = "B-1a","B1b"= "B-1b"))
write_csv(bcr.light, "../bcr.light.csv")



############## ############## ############## ############## ############## ############## ############## ############## ############## 
############## ############## ############## ##############  ISOTYPES AND UNIQUE CDR3 ############## ############## ############## 
############## ############## ############## ############## ############## ############## ############## ############## ############## 
bcr.light$Class <- ifelse(str_detect(bcr.light$C, "IGHG"), "IGHG", bcr.light$C)
st.Ig <- bcr.light%>% 
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
all_Bs <- st.Ig  %>% dplyr::filter(Class != ".", cell_type != "B2") %>%
  ggplot(aes(x = condition, y = Num.Ig, fill = condition))+ 
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  labs(y = "Normalized Ig Abundance",fill = "Ig")+
  scale_fill_manual(values =  cls)+
  theme_prism(base_size = 10) +facet_wrap(~Class + cell_type, scale = "free_y") +
  scale_y_continuous(expand = expansion(mult = 0.3)) +
  stat_compare_means(method = "wilcox.test", hide.ns = T, vjust = -.7, label.x.npc  = 0.5, label = "p.signif", size = 5, fontface = "bold") + geom_jitter(alpha = 0.5) 

IGKC <- st.Ig %>% mutate(cell_type = recode_factor(cell_type, "B1a" = "B-1a","B-1b"= "B-1b")) %>% 
  dplyr::filter(Class == "IGKC", cell_type != "B2") %>%
  ggplot(aes(x = condition, y = Num.Ig, fill = condition))+ 
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  labs(y = "Normalized Ig Abundance",fill = "Ig")+
  scale_fill_manual(values =  cls)+
  theme_prism(base_size = 10) +facet_wrap(~cell_type) +
  scale_y_continuous(expand = expansion(mult = 0.3)) + ggtitle("IGKC") +
  stat_compare_means(method = "wilcox.test", hide.ns = T, vjust = -.7, label.x.npc  = 0.5, label = "p.signif", size = 5, fontface = "bold") + geom_jitter(alpha = 0.5) 
IGLC <- st.Ig%>% mutate(cell_type = recode_factor(cell_type, "B1a" = "B-1a","B1b"= "B-1b")) %>% 
  dplyr::filter(Class == "IGLC", cell_type != "B2") %>%
  ggplot(aes(x = condition, y = Num.Ig, fill = condition))+ 
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  labs(y = "Normalized Ig Abundance",fill = "Ig")+
  scale_fill_manual(values =  cls)+
  theme_prism(base_size = 10) +facet_wrap(~cell_type) +
  scale_y_continuous(expand = expansion(mult = 0.3)) + ggtitle("IGLC") +
  stat_compare_means(method = "wilcox.test", hide.ns = T, vjust = -.7, label.x.npc  = 0.5, label = "p.signif", size = 5, fontface = "bold") + geom_jitter(alpha = 0.5) 
plot <- ggarrange(IGKC, IGLC, nrow = 2, common.legend = TRUE, legend = "right")+bgcolor("white") +border("white")
ggsave( plot = plot,
        filename = file_save("all_igs_trust4.tiff"), 
        height = 5, width = 4.12, units = "in", dpi = 300)

all_Bs
ggsave( plot = all_Bs,
        filename = file_save("all_igs_trust4.tiff"), 
        height = 4, width = 4.12, units = "in", dpi = 300)
all_Bs2 <- st.Ig %>% mutate(Class = factor(Class, levels = c("IGKC", "IGLC"))) %>% 
  dplyr::filter(Class != ".", cell_type != "B2") %>%
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
        filename = file_save("all_igs_trust4.3.tiff"), 
        height = 4.6, width = 4.12, units = "in", dpi = 300)

B1a.bcr.light <- bcr.light %>% filter(cell_type == "B-1a")
dat <- aggregate(CDR3aa ~ sample + condition, B1a.bcr.light, function(x) length(unique(x))) 
B1a_gc <- CompareGroups(dat,"CDR3aa","Number of Unique CDR3", "B-1a")

B1b.bcr.light <- bcr.light %>% filter(cell_type == "B-1b")
dat <- aggregate(CDR3aa ~ sample + condition, B1b.bcr.light, function(x) length(unique(x))) 
B1b_gc <- CompareGroups(dat,"CDR3aa","Number of Unique CDR3", "B-1b")
B1b_gc <- B1b_gc + ylab("")

B2.bcr.light <- bcr.light %>% filter(cell_type == "B2")
dat <- aggregate(CDR3aa ~ sample + condition, B2.bcr.light, function(x) length(unique(x))) 
B2_gc <- CompareGroups(dat,"CDR3aa","Number of Unique CDR3", "B2")


pl <- ggarrange(B1a_gc, B1b_gc,  nrow = 1, common.legend = TRUE, legend = "bottom")
pl
ggsave(plot = pl, file = file_save("unique_CDR3.pdf"),h =2.54, w = 4.54, dpi = 300)
.V <- bcr.light %>% dplyr::group_by(sample, condition, cell_type, V) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
.D <- bcr.light %>% dplyr::group_by(sample, condition, cell_type, D) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
.J <- bcr.light %>% dplyr::group_by(sample, condition, cell_type, J) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
write_csv(.V, "V.csv")
write_csv(.D, "D.csv")
write_csv(.J, "J.csv")
write_csv(bcr.light, "heavy_chain.csv")

############## ############## ############## ############## ############## ############## ############## ############## 
############## CDR3
############## ############## ############## ############## ############## ############## ############## ############## 

cdr3 <- bcr.light %>% dplyr::group_by(sample, condition, cell_type, CDR3aa) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
tmp1 <- cdr3 %>% complete(nesting(smp, condition, cell_type), CDR3aa) %>% replace(is.na(.), 0)
B1a_p_val_cdr3 <- tmp1 %>% filter(cell_type == "B-1a") %>% 
  dplyr::group_by(CDR3aa) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(CDR3aa, Wilcox = w$p.value)
sigs <- B1a_p_val_cdr3 %>% filter(Wilcox < 0.05, CDR3aa != "out_of_frame") %>% pull(CDR3aa)
cdr3_in_B1a <- tmp1 %>% filter(cell_type == "B-1a", CDR3aa != "out_of_frame") %>% group_by(condition, CDR3aa) %>% 
  summarize(avg_freq = mean(freq)) %>% ungroup %>% 
  arrange(desc(avg_freq)) %>% select(CDR3aa) %>%
  unique %>% slice(1:25) %>%  pull(CDR3aa)
tmp1 %>% filter(cell_type == "B-1a", CDR3aa != "out_of_frame") %>% 
  arrange(desc(freq)) %>% pull(CDR3aa) %>% unique()

B1a_CDR3aa <- tmp1 %>% 
  filter(cell_type == "B-1a", CDR3aa %in% cdr3_in_B1a,  CDR3aa != "out_of_frame") %>% 
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
       filename = file_save(
                         "B1a_CDR3aa_trust4.tiff"),
       height = 5.84, width = 6.84, unit = "in", dpi =300)

cdr3 <- bcr.light %>% dplyr::group_by(sample, condition, cell_type, CDR3aa) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
tmp1 <- cdr3 %>% complete(nesting(smp, condition, cell_type), CDR3aa) %>% replace(is.na(.), 0)
B1b_p_val_cdr3 <- tmp1 %>% filter(cell_type == "B-1b") %>% 
  dplyr::group_by(CDR3aa) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(CDR3aa, Wilcox = w$p.value)
sigs <- B1b_p_val_cdr3 %>% filter(Wilcox < 0.05, CDR3aa != "out_of_frame") %>% pull(CDR3aa)
cdr3_in_B1b <- tmp1 %>% filter(cell_type == "B-1b", CDR3aa != "out_of_frame") %>% group_by(condition, CDR3aa) %>% 
  summarize(avg_freq = mean(freq)) %>% ungroup %>% 
  arrange(desc(avg_freq)) %>% select(CDR3aa) %>%
  unique %>% slice(1:25) %>%  pull(CDR3aa)
tmp1 %>% filter(cell_type == "B-1b", CDR3aa != "out_of_frame") %>% 
  arrange(desc(freq)) %>% pull(CDR3aa) %>% unique()

B1b_CDR3aa <- tmp1 %>% 
  filter(cell_type == "B-1b", CDR3aa %in% cdr3_in_B1b,  CDR3aa != "out_of_frame") %>% 
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
       filename = file_save(
                         "B1b_CDR3aa_trust4.tiff"),
       height = 5.84, width = 6.84, unit = "in", dpi =300)
###############################
pl1 <- bcr.light %>% 
  filter(cell_type == "B-1a",CDR3aa != "out_of_frame", condition == "WT")  %>%  slice(1:20) %>% 
  ggplot(aes(reorder(CDR3aa, frequency), frequency, fill = frequency)) + 
  geom_bar(stat = "identity") +
  coord_flip() + 
  scale_fill_gradient(high="orange", low="blue") +  coord_flip() +
  scale_y_continuous(labels = scales::percent, limits = c(0, .4)) + 
  xlab("CDR3aa")  + 
  ylab("% of all CDR3") + theme_prism(base_size = 8)

pl2 <- bcr.light %>% 
  filter(cell_type == "B-1a",CDR3aa != "out_of_frame", condition == "KO")  %>%  slice(1:20) %>% 
  ggplot(aes(reorder(CDR3aa, frequency), frequency, fill = frequency)) + 
  geom_bar(stat = "identity") +
  coord_flip() + 
  scale_fill_gradient(high="orange", low="blue") +  coord_flip() +
  scale_y_continuous(labels = scales::percent, limits = c(0, .4)) + 
  xlab("CDR3aa")  + 
  ylab("% of all CDR3") + theme_prism(base_size = 8)

plot <- ggarrange(pl1, pl2, nrow = 1, common.legend = TRUE, legend = "right")+bgcolor("white") +border("white")
ggsave(file_save("CDR3aa_B1a.png"), plot = plot,  height = 4.5, width = 8, unit = "in")

pl1 <- bcr.light %>% 
  filter(cell_type == "B-1b",CDR3aa != "out_of_frame", condition == "WT")  %>%  slice(1:20) %>% 
  ggplot(aes(reorder(CDR3aa, frequency), frequency, fill = frequency)) + 
  geom_bar(stat = "identity") +
  coord_flip() + 
  scale_fill_gradient(high="orange", low="blue") +  coord_flip() +
  scale_y_continuous(labels = scales::percent, limits = c(0, .4)) + 
  xlab("CDR3aa")  + 
  ylab("% of all CDR3") + theme_prism(base_size = 8)

pl2 <- bcr.light %>% 
  filter(cell_type == "B-1b",CDR3aa != "out_of_frame", condition == "KO")  %>%  slice(1:20) %>% 
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
############## ############## ##############  IGKCs
############## ############## ############## ############## ############## ############## ############## ############## ############## 
IGKCs <- bcr.light %>% filter(cell_type == "B-1a", Class == "IGKC")
IGKCs <- IGKCs %>% dplyr::group_by(sample, condition, cell_type, CDR3aa) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
tmp1 <- IGKCs %>% complete(nesting(smp, condition, cell_type), CDR3aa) %>% replace(is.na(.), 0)
IGKCs_p_val_cdr3 <- tmp1 %>% filter(cell_type == "B-1a") %>% 
  dplyr::group_by(CDR3aa) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(CDR3aa, Wilcox = w$p.value)
sigs <- IGKCs_p_val_cdr3 %>% filter(Wilcox < 0.05, CDR3aa != "out_of_frame") %>% pull(CDR3aa)
IGKCs_in_B1b <- tmp1 %>% filter(cell_type == "B-1a", CDR3aa != "out_of_frame") %>% group_by(condition, CDR3aa) %>% 
  summarize(avg_freq = mean(freq)) %>% ungroup %>% 
  arrange(desc(avg_freq)) %>% select(CDR3aa) %>%
  unique %>% slice(1:25) %>%  pull(CDR3aa)
tmp1 %>% filter(cell_type == "B-1a", CDR3aa != "out_of_frame") %>% 
  arrange(desc(freq)) %>% pull(CDR3aa) %>% unique()


IGKCs <- tmp1 %>% 
  filter(cell_type == "B-1a", CDR3aa %in% IGKCs_in_B1b,  CDR3aa != "out_of_frame") %>% 
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
IGKCs



###### ALternative V analysis
all_V <- bcr.light
all_V$V <- gsub("\\-.*", "", bcr.light$V)
all_V$V <- gsub("\\*.*", "", all_V$V)

all_V$V <- gsub("S.*", "", all_V$V)
all_V$V <- gsub("IGKV", "VK", all_V$V)
all_V$V <- gsub("IGLV", "VL", all_V$V)

v <- all_V %>% filter(V != ".") %>% dplyr::group_by(sample, condition, cell_type, V) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct

tmp1 <- v %>% complete(nesting(smp, condition, cell_type), V) %>% replace(is.na(.), 0)
B1a_p_val <- tmp1 %>% filter(cell_type == "B-1a") %>% 
  dplyr::group_by(V) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(V, Wilcox = w$p.value)

sigs <- B1a_p_val %>% filter(Wilcox < 0.05) %>% pull(V)
B1a_V <- tmp1 %>% filter(cell_type == "B-1a") %>% ggplot(aes(x = reorder(V, -freq), y = freq, fill = condition )) +
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
ggsave(plot = B1a_V, filename = file_save("B1a_v_VH.tiff"),
       height = 3.3, width = 6.6, unit = "in", dpi =300)



tmp1 <- v %>% complete(nesting(smp, condition, cell_type), V) %>% replace(is.na(.), 0)
B1b_p_val <- tmp1 %>% filter(cell_type == "B-1b") %>% 
  dplyr::group_by(V) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(V, Wilcox = w$p.value)

sigs <- B1b_p_val %>% filter(Wilcox < 0.05) %>% pull(V)
B1b_V <- tmp1 %>% filter(cell_type == "B-1b") %>% ggplot(aes(x = reorder(V, -freq), y = freq, fill = condition )) +
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
ggsave(plot = B1b_V, filename =file_save("B1b_v_VH.tiff"),
       height = 3.3, width = 6.6, unit = "in", dpi =300)
######
########## J
all_J <- bcr.light
all_J$J <- gsub("\\*.*", "", bcr.light$J)
all_J$J <- gsub("IGLJ", "JL", all_J$J)
all_J$J <- gsub("IGKJ", "JK", all_J$J)

J <- all_J %>% filter(J != ".") %>% dplyr::group_by(sample, condition, cell_type, J) %>%
  summarize(freq = sum(frequency)) %>%
  ungroup %>%  mutate(smp = sample) %>% select(-sample) %>% distinct
tmp1 <- J %>% complete(nesting(smp, condition, cell_type), J) %>% replace(is.na(.), 0)
B1a_p_val <- tmp1 %>% filter(cell_type == "B-1a") %>% 
  dplyr::group_by(J) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(J, Wilcox = w$p.value)

sigs <- B1a_p_val %>% filter(Wilcox < 0.05) %>% pull(J)
B1a_J <- tmp1 %>% filter(cell_type == "B-1a") %>% ggplot(aes(x = J, y = freq, fill = condition )) +
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
ggsave(plot = B1a_J, filename = file_save("B1a_J_JH.tiff"),
       height = 3.3, width = 6.6, unit = "in", dpi =300)


B1b_p_val <- tmp1 %>% filter(cell_type == "B-1b") %>% 
  dplyr::group_by(J) %>% 
  do(w = wilcox.test(freq ~ condition, data = ., exact = F, paired = F)) %>%      
  summarise(J, Wilcox = w$p.value)

sigs <- B1b_p_val %>% filter(Wilcox < 0.05) %>% pull(J)
B1b_J <- tmp1 %>% filter(cell_type == "B-1b") %>% ggplot(aes(x = J, y = freq, fill = condition )) +
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
ggsave(plot = B1b_J, filename = file_save("B1b_J_JH.tiff"),
       height = 3.3, width = 6.6, unit = "in", dpi =300)




###########
########### V presence in replicated sequences
colo <- c('#4363d8', '#f58231','#e6194b', '#3cb44b', '#ffe119',  '#911eb4',
          '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
          '#9a6324', '#fffac8', '#800000', '#aaffc3',
          '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

bcr.light.v <- bcr.light

CDR3aa <- bcr.light.v %>% filter(cell_type == "B-1a", CDR3aa != "out_of_frame") %>%
  select(cell_type, condition,CDR3aa, CDR3nt, count) %>% 
  group_by(cell_type, condition) %>%
  mutate(seq_n = sum(count)) %>%
  ungroup %>% 
  group_by(cell_type, condition, CDR3aa)  %>% 
  summarize(c_CDR3 = sum(count), f_CDR3 = sum(count)/seq_n ) %>%
  ungroup %>% distinct %>% arrange(-f_CDR3) %>% unique

CDR3aa.V <- bcr.light.v %>% filter(cell_type == "B-1a",  CDR3aa != "out_of_frame") %>%
  select(cell_type, condition,CDR3aa, CDR3nt, count, V) %>% 
  group_by(cell_type, condition) %>%
  mutate(seq_n = sum(count)) %>%
  ungroup %>% 
  group_by(cell_type, condition, CDR3aa, V)  %>% 
  summarize(c_CDR3 = sum(count), f_CDR3 = sum(count)/seq_n ) %>%
  ungroup %>% distinct %>% arrange(-f_CDR3) %>% ungroup %>% unique
  

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
ggsave(plot = cd1, filename = file_save("piechart_B1a_WT_cdr3.tiff"), dpi = 300, w = 7, h = 4.7)

cdr3_props_second_plot <- CDR3_props %>% 
  filter(condition == "KO") %>% arrange(c)

cdr3_props_second_plot$CDR3 <- fct_reorder(cdr3_props_second_plot$CDR3, -(cdr3_props_second_plot$c) )

cdr3_props_second_plot$CDR3 <- relevel(cdr3_props_second_plot$CDR3, ref = "not replicated")

cd2 <- cdr3_props_second_plot %>% filter(condition == "KO") %>% ggplot(aes(x = "", y=c, fill = CDR3)) +
  geom_bar(stat="identity", width=0.001, color="white") + labs(fill = "CDR3 AA sequence (B-1a KO)")+
  coord_polar("y", start=0) +
  scale_fill_manual(values = colo) +
  theme_classic(base_size = 14) %+replace% theme(line = element_blank(), axis.text = element_blank()) +
  ylab("") + xlab("") 
cd2
ggsave(plot = cd2, filename = file_save("piechart_B1a_KO_cdr3.tiff"), dpi = 300, w = 7, h = 4.7)
CDR3_props <- CDR3aa.V %>% select(CDR3aa, condition, f_CDR3, c_CDR3) %>% 
  group_by(condition) %>%  unique %>% ungroup
CDR3_props$rep<- ifelse(CDR3_props$f_CDR3 > 0.01, "replicated", "not replicated")
table(CDR3_props$condition, CDR3_props$rep)
f_test <- CDR3_props %>% mutate(f = round(f_CDR3 * 100))
cond1 <- rep(f_test$condition ,f_test$c_CDR3)
cond2 <- rep(f_test$rep ,f_test$c_CDR3)
cond <- tibble(condition = cond1, rep = cond2)
table(cond$condition, cond$rep)
fisher.test(cond$condition, cond$rep)
table(CDR3_props$condition, CDR3_props$rep)
c_test <- CDR3_props %>% mutate(c = round(c_CDR3))

cond1 <- rep(c_test$condition ,c_test$c)
cond2 <- rep(c_test$rep ,c_test$c)
cond <- tibble(condition = cond1, rep = cond2)
table(cond$condition, cond$rep)
rep_test <-CDR3_props  %>% group_by(condition) %>% mutate(total = n()) %>% ungroup %>% 
  group_by(condition, rep, total) %>%  summarize(n = n()) %>% mutate(perc = n/total)

rep_test_all <-CDR3_props  %>% group_by(condition) %>% mutate(scaling = sum(c_CDR3)/1000000) %>% 
  mutate(norm_counts = c_CDR3/scaling) %>% ungroup %>% 
  group_by(condition, rep) %>%
  summarize(n = round(sum(norm_counts))) %>% unique

test <- matrix(c(rep_test_all$n[1],rep_test_all$n[3],rep_test_all$n[2], rep_test_all$n[4]), nrow = 2)
chisq.test(test)
chisq.test( table(CDR3_props$condition, CDR3_props$rep))
table(CDR3_props$condition, CDR3_props$rep)
test

CDR3aa <- bcr.light.v %>% filter(cell_type == "B-1b", CDR3aa != "out_of_frame") %>%
  select(cell_type, condition,CDR3aa, CDR3nt, count) %>% 
  group_by(cell_type, condition) %>%
  mutate(seq_n = sum(count)) %>%
  ungroup %>% 
  group_by(cell_type, condition, CDR3aa)  %>% 
  summarize(c_CDR3 = sum(count), f_CDR3 = sum(count)/seq_n ) %>%
  ungroup %>% distinct %>% arrange(-f_CDR3) %>% unique

CDR3aa.V <- bcr.light.v %>% filter(cell_type == "B-1b", CDR3aa != "out_of_frame") %>%
  select(cell_type, condition,CDR3aa, CDR3nt, count, V) %>% 
  group_by(cell_type, condition) %>%
  mutate(seq_n = sum(count)) %>%
  ungroup %>% 
  group_by(cell_type, condition, CDR3aa, V)  %>% 
  summarize(c_CDR3 = sum(count), f_CDR3 = sum(count)/seq_n ) %>%
  ungroup %>% distinct %>% arrange(-f_CDR3) %>% ungroup %>% unique


CDR3_props <- CDR3aa.V %>% select(CDR3aa, condition, f_CDR3, c_CDR3) %>% 
  group_by(condition) %>%  unique %>% ungroup
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
ggsave(plot = cd1, filename = file_save("piechart_B1b_WT_cdr3.tiff"), dpi = 300, w = 7, h = 4.57)
cd2 <- CDR3_props %>% filter(condition == "KO") %>% ggplot(aes(x = "", y=c, fill = reorder(CDR3, -c))) +
  geom_bar(stat="identity", width=0.001, color="white") + labs(fill = "CDR3 AA sequence (B-1b KO)")+
  coord_polar("y", start=0) +
  scale_fill_manual(values = colo) +
  theme_classic(base_size = 14) %+replace% theme(line = element_blank(), axis.text = element_blank()) +
  ylab("") + xlab("") 
ggsave(plot = cd2, filename = file_save("piechart_B1b_KO_cdr3.tiff"), dpi = 300, w = 7, h = 5.2)

CDR3_props <- CDR3aa.V %>% select(CDR3aa, condition, f_CDR3, c_CDR3) %>% 
  group_by(condition) %>%  unique %>% ungroup
CDR3_props$rep<- ifelse(CDR3_props$f_CDR3 > 0.01, "replicated", "not replicated")
table(CDR3_props$condition, CDR3_props$rep)
f_test <- CDR3_props %>% mutate(f = round(f_CDR3 * 100))
cond1 <- rep(f_test$condition ,f_test$c_CDR3)
cond2 <- rep(f_test$rep ,f_test$c_CDR3)
cond <- tibble(condition = cond1, rep = cond2)
table(cond$condition, cond$rep)
fisher.test(cond$condition, cond$rep)
table(CDR3_props$condition, CDR3_props$rep)
c_test <- CDR3_props %>% mutate(c = round(c_CDR3))

cond1 <- rep(c_test$condition ,c_test$c)
cond2 <- rep(c_test$rep ,c_test$c)
cond <- tibble(condition = cond1, rep = cond2)
table(cond$condition, cond$rep)
rep_test <-CDR3_props  %>% group_by(condition) %>% mutate(total = n()) %>% ungroup %>% 
  group_by(condition, rep, total) %>%  summarize(n = n()) %>% mutate(perc = n/total)

rep_test_all <-CDR3_props  %>% group_by(condition) %>% mutate(scaling = sum(c_CDR3)/1000000) %>% 
  mutate(norm_counts = c_CDR3/scaling) %>% ungroup %>% 
  group_by(condition, rep) %>%
  summarize(n = round(sum(norm_counts))) %>% unique

test <- matrix(c(rep_test_all$n[1],rep_test_all$n[3],rep_test_all$n[2], rep_test_all$n[4]), nrow = 2)
chisq.test(test)
chisq.test( table(CDR3_props$condition, CDR3_props$rep))
table(CDR3_props$condition, CDR3_props$rep)
test
