library(tidyverse)
suppressMessages(library(ggpubr))
library(ggpubr)
#library(rstatix)
library(tidyr)
library(ggprism)
file_save <- function(.fname){
  if (!dir.exists("./plots/")){
    dir.create("./plots/")
  }else{
    print("dir exists")
  }
  filename = paste0("./plots/", .fname)
  return(filename)
}
samples <- read_csv("samples.csv")
rownames(samples) <- samples$library
bcr.heavy <- read_csv("bcr_heavy.csv.gz")
bcr.heavy$condition <- factor(bcr.heavy$condition, levels = c("WT", "KO"))
bcr.heavy <- bcr.heavy %>% mutate(cell_type = recode_factor(cell_type, "B1a" = "B-1a","B1b"= "B-1b"))

colo <- c('#4363d8', '#f58231','#e6194b', '#3cb44b', '#ffe119',  '#911eb4',
          '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
          '#9a6324', '#fffac8', '#800000', '#aaffc3',
          '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
bcr.heavy.v <- bcr.heavy
bcr.heavy.v$V <- gsub("\\-.*", "", bcr.heavy.v$V)
bcr.heavy.v$V <- gsub("S.*", "", bcr.heavy.v$V)
bcr.heavy.v$V <- gsub("IGHV", "VH", bcr.heavy.v$V)
bcr.heavy.v$J<- gsub("\\*..", "", bcr.heavy.v$J)
bcr.heavy.v$J <- gsub("IGHJ", "JH", bcr.heavy.v$J)
CDR3aa <- bcr.heavy.v %>% filter(cell_type == "B-1a") %>%
  select(cell_type, condition,CDR3aa, CDR3nt, count, V, J) %>% 
  group_by(cell_type, condition) %>%
  mutate(seq_n = sum(count)) %>%
  ungroup %>% 
  group_by(cell_type, condition, CDR3aa, V, J)  %>% 
  dplyr::summarize(c_CDR3 = sum(count), f_CDR3 = sum(count)/seq_n ) %>%
  ungroup %>% distinct %>% arrange(-f_CDR3)




CDR3_props <- CDR3aa %>% select(CDR3aa, condition, f_CDR3, c_CDR3, V, J) %>% 
  group_by(condition) %>%  unique %>% ungroup
table(CDR3_props$condition, CDR3_props$f_CDR3 > 0.01)
seq_clones_WT <- sum((CDR3_props[CDR3_props$f_CDR3 > 0.01 & CDR3_props$condition == "WT",])$f_CDR3)
seq_clones_KO <- sum((CDR3_props[CDR3_props$f_CDR3 > 0.01 & CDR3_props$condition == "KO",])$f_CDR3)
CDR3_props$CDR3 <- ifelse(CDR3_props$f_CDR3 > 0.01, CDR3_props$CDR3aa, "not replicated")
CDR3_props <- CDR3_props %>%  group_by(condition, CDR3, V, J) %>%  
  dplyr::summarize(f = sum(f_CDR3), c = sum(c_CDR3))

CDR3_props$CDR3_full <- ifelse(CDR3_props$CDR3 == "not replicated", "not replicated", 
       paste0(CDR3_props$CDR3, " (", CDR3_props$V, " + ", CDR3_props$J, ")"))
CDR3_props <- CDR3_props %>% select(condition, f, c, CDR3_full) %>% ungroup
CDR3_props <- CDR3_props %>%  group_by(condition, CDR3_full) %>%  
  dplyr::summarize(f = sum(f), c = sum(c)) %>% ungroup 


cd1 <- CDR3_props %>% filter(condition == "WT") %>% 
  ggplot(aes(x = "", y=c,fill = reorder(CDR3_full, -c))) +
  geom_bar(stat="identity", width=0.001, color="white") + labs(fill = "CDR3 AA sequence (B-1a WT)") +
  coord_polar("y", start=0) + 
  scale_fill_manual(values = colo) +
  theme_classic(base_size = 14) %+replace% theme(line = element_blank(), axis.text = element_blank()) +
  ylab("") + xlab("") 

ggsave(plot = cd1, filename = file_save("piechart_B1a_WT_cdr3_VJ.pdf"), dpi = 300, w = 7, h = 4.5)

colo <- c('#4363d8', '#3cb44b', '#ffe119','#e6194b', '#f58231',  '#911eb4',
          '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
          '#9a6324', '#fffac8', '#800000', '#aaffc3',
          '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
cd2 <- CDR3_props %>% filter(condition == "KO") %>% ggplot(aes(x = "", y=c, fill = reorder(CDR3_full, -c))) +
  geom_bar(stat="identity", width=0.001, color="white") + labs(fill = "CDR3 AA sequence (B-1a KO)")+
  coord_polar("y", start=0) +
  scale_fill_manual(values = colo) +
  theme_classic(base_size = 14) %+replace% theme(line = element_blank(), axis.text = element_blank()) +
  ylab("") + xlab("") 
ggsave(plot = cd2, filename = file_save("piechart_B1a_KO_cdr3_VJ.pdf"), dpi = 300, w = 7, h = 4.5)


colo <- c('#4363d8', '#f58231','#e6194b', '#3cb44b', '#ffe119',  '#911eb4',
          '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
          '#9a6324', '#fffac8', '#800000', '#aaffc3',
          '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

#### same for B1b
CDR3aa <- bcr.heavy.v %>% filter(cell_type == "B-1b") %>%
  select(cell_type, condition,CDR3aa, CDR3nt, count, V, J) %>% 
  group_by(cell_type, condition) %>%
  mutate(seq_n = sum(count)) %>%
  ungroup %>% 
  group_by(cell_type, condition, CDR3aa, V, J)  %>% 
  dplyr::summarize(c_CDR3 = sum(count), f_CDR3 = sum(count)/seq_n ) %>%
  ungroup %>% distinct %>% arrange(-f_CDR3)




CDR3_props <- CDR3aa %>% select(CDR3aa, condition, f_CDR3, c_CDR3, V, J) %>% 
  group_by(condition) %>%  unique %>% ungroup %>% mutate(f_CDR3 = round(f_CDR3, digits = 4))
table(CDR3_props$condition, CDR3_props$f_CDR3 >= 0.01)
seq_clones_WT <- sum((CDR3_props[CDR3_props$f_CDR3 >= 0.01 & CDR3_props$condition == "WT",])$f_CDR3)
seq_clones_KO <- sum((CDR3_props[CDR3_props$f_CDR3 >= 0.01 & CDR3_props$condition == "KO",])$f_CDR3)
CDR3_props$CDR3 <- ifelse(CDR3_props$f_CDR3 >= 0.01, CDR3_props$CDR3aa, "not replicated")
CDR3_props <- CDR3_props %>%  group_by(condition, CDR3, V, J) %>%  
  dplyr::summarize(f = sum(f_CDR3), c = sum(c_CDR3))

CDR3_props$CDR3_full <- ifelse(CDR3_props$CDR3 == "not replicated", "not replicated", 
                               paste0(CDR3_props$CDR3, " (", CDR3_props$V, " + ", CDR3_props$J, ")"))
CDR3_props <- CDR3_props %>% select(condition, f, c, CDR3_full) %>% ungroup
CDR3_props <- CDR3_props %>%  group_by(condition, CDR3_full) %>%  
  dplyr::summarize(f = sum(f), c = sum(c)) %>% ungroup 


cd1 <- CDR3_props %>% filter(condition == "WT") %>% 
  ggplot(aes(x = "", y=c,fill = reorder(CDR3_full, -c))) +
  geom_bar(stat="identity", width=0.001, color="white") + labs(fill = "CDR3 AA sequence (B-1b WT)") +
  coord_polar("y", start=0) + 
  scale_fill_manual(values = colo) +
  theme_classic(base_size = 14) %+replace% theme(line = element_blank(), axis.text = element_blank()) +
  ylab("") + xlab("") 

ggsave(plot = cd1, filename = file_save("piechart_B1b_WT_cdr3_VJ.pdf"), dpi = 300, w = 7, h = 4.5)

cd2 <- CDR3_props %>% filter(condition == "KO") %>% ggplot(aes(x = "", y=c, fill = reorder(CDR3_full, -c))) +
  geom_bar(stat="identity", width=0.001, color="white") + labs(fill = "CDR3 AA sequence (B-1b KO)")+
  coord_polar("y", start=0) +
  scale_fill_manual(values = colo) +
  theme_classic(base_size = 14) %+replace% theme(line = element_blank(), axis.text = element_blank()) +
  ylab("") + xlab("") 
ggsave(plot = cd2, filename = file_save("piechart_B1b_KO_cdr3_VJ.pdf"), dpi = 300, w = 7, h = 4.5)

