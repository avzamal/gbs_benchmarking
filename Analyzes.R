library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(car)
library(rstatix)
library(emmeans)
library(agricolae)
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path) 
setwd(main_dir)
options(scipen=999)

# Load data and preparation
NRD <- read.csv("benchmarking_table_NRD_str.txt", sep="")
r_square <- read.csv("benchmarking_table_r_square_str.txt", sep="")
tstv <- read.csv("benchmarking_table_tstv_str.txt", sep="")

full_pairs <- left_join(left_join(NRD, r_square, by="pair"),tstv,by="pair")
full_pairs <- separate(full_pairs,pair,sep = "_",into = c("Enzyme_1","Aligner_1","Caller_1","Enzyme_2","Aligner_2","Caller_2")) %>% mutate(Aligner_1 = if_else(Aligner_1 == "str", "Strobealign", Aligner_1)) %>% 
  mutate(Aligner_2 = if_else(Aligner_2 == "str", "Strobealign", Aligner_2)) %>%
  unite("Label_1",1:3,sep = "_",remove=F) %>% 
  unite("Label_2",5:7,sep = "_",remove=F) %>%
  unite("Combination_2",7:8,sep = "_",remove=F) %>%
  unite("Combination_1",3:4,sep = "_",remove=F) %>% 
  unite("Enzyme_pair",c(2,7),sep = "_",remove=F) 
View(full_pairs)

calling <- read.table("benchmarking_table_calling_str.txt", quote="\"", comment.char="")
WGS_GBS <- read.table("benchmarking_table_WGS_GBS_str.txt", quote="\"", comment.char="")
not_in_genes <- read.table("benchmarking_table_not_in_genes_str.txt", quote="\"", comment.char="")
not_in_genes <- separate(not_in_genes, col = 2,into = c("Label","Rm"),sep = "_not") %>% select(V1,Label)
colnames(not_in_genes)[1] <- "Not_in_genes"
not_in_repeats <- read.table("benchmarking_table_not_in_repeats_str.txt", quote="\"", comment.char="")
not_in_repeats <- separate(not_in_repeats, col = 2,into = c("Label","Rm"),sep = "_not") %>% select(V1,Label)
colnames(not_in_repeats)[1] <- "Not_in_repeats"

individual_stats <- cbind(calling,WGS_GBS)
colnames(individual_stats) <- c("total_snps","file","WGS_intersection")
individual_stats <- separate(individual_stats,file,sep = "_",into = c("Enzyme","Aligner","Caller","ext")) %>%
  select(-ext) %>% 
  unite("Label",2:4,sep = "_",remove=F) %>% 
  unite("Combination",4:5,sep = "_",remove=F) %>% 
  select(Label,Combination,Enzyme,Aligner,Caller,total_snps,WGS_intersection) %>% 
  mutate(ratio = WGS_intersection/total_snps) %>% 
  left_join(not_in_genes) %>% 
  mutate(in_genes = total_snps - Not_in_genes) %>% 
  mutate(ratio_genes = in_genes/total_snps) %>% 
  left_join(not_in_repeats) %>% 
  mutate(in_repeats = total_snps - Not_in_repeats) %>% 
  mutate(ratio_repeats = in_repeats/total_snps) %>% 
  mutate(Aligner = if_else(Aligner == "str", "Strobealign", Aligner)) %>% 
  unite("Label",2:4,sep = "_",remove=F) %>% 
  unite("Combination",4:5,sep = "_",remove=F)

individual_stats$Enzyme <- factor(individual_stats$Enzyme,labels = c("ApeKI","HindIII-NlaIII","PstI-MspI"))
individual_stats <- individual_stats %>% group_by(Enzyme) %>% mutate(SNPs_ratio = total_snps/max(total_snps)) %>% ungroup()
individual_stats %>% group_by(Caller) %>% summarise(mean(ratio))
individual_stats %>% group_by(Caller) %>% summarise(mean(SNPs_ratio))
View(individual_stats)

# Preparation data for plotting
# Subset of pairs which came from the same restriction enzyme
full_pairs_per_enzyme <- filter(full_pairs,Enzyme_1 == Enzyme_2)

# Plots of aligner with different enzymes
(a5 <- ggplot(filter(individual_stats), aes(Enzyme,total_snps,fill=Aligner))+
    geom_boxplot() +
    labs(x = "Enzyme set",
         y="Total number of SNPs")+
    theme_bw(base_size = 16))
ggsave("Aligner_Total.png",a5)
(d5 <- ggplot(filter(individual_stats), aes(Enzyme, ratio, fill=Aligner))+
    geom_boxplot() +
    labs(x = "Enzyme set",
         y="SNPs presented in WGS")+
    theme_bw(base_size = 16))
ggsave("Aligner_Ratio.png",d5)

# Plots without comparison ApeKI Aligners
(a <- ggplot(filter(individual_stats, Enzyme == "ApeKI"), aes(Aligner,total_snps))+
  geom_boxplot() +
  labs(x = "Aligner",
       y="Total number of SNPs after filtering")+
  theme_bw(base_size = 16))
ggsave("ApeKI_Aligner_Total.png",a)


(c <- ggplot(filter(individual_stats, Enzyme == "ApeKI"&total_snps<300000), aes(Aligner, WGS_intersection))+
    geom_boxplot() +
    labs(x = "Aligner",
         y="SNPs presented in WGS")+
    theme_bw(base_size = 16))
ggsave("ApeKI_Aligner_WGS.png",c)

(d <- ggplot(filter(individual_stats, Enzyme == "ApeKI"&total_snps<300000), aes(Aligner, ratio))+
    geom_boxplot() +
    labs(x = "Aligner",
         y="SNPs presented in WGS")+
    theme_bw(base_size = 16))
ggsave("ApeKI_Aligner_Ratio.png",d)

# Plots without comparison ApeKI callers
(e <- ggplot(filter(individual_stats, Enzyme == "ApeKI"), aes(Caller,total_snps))+
    geom_boxplot() +
    labs(x = "Caller",
         y="Total number of SNPs after filtering")+
    theme_bw(base_size = 16)+
    theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust=1)))
ggsave("ApeKI_Caller_Total.png",e)

(g <- ggplot(filter(individual_stats, Enzyme == "ApeKI"&total_snps<300000), aes(Caller, WGS_intersection))+
    geom_boxplot() +
    labs(x = "Aligner",
         y="SNPs presented in WGS")+
    theme_bw(base_size = 16)+
    theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust=1)))
ggsave("ApeKI_Caller_WGS.png",g)

(h <- ggplot(filter(individual_stats, Enzyme == "ApeKI"&total_snps<300000), aes(Caller, ratio))+
    geom_boxplot() +
    labs(x = "Caller",
         y="SNPs presented in WGS")+
    theme_bw(base_size = 16)+
    theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust=1)))
ggsave("ApeKI_Caller_Ratio.png",h)

# Plots without comparison PstI Aligners
(a1 <- ggplot(filter(individual_stats, Enzyme == "PstI-MspI"), aes(Aligner,total_snps))+
    geom_boxplot() +
    labs(x = "Aligner",
         y="Total number of SNPs after filtering")+
    theme_bw(base_size = 16))
ggsave("PstI_Aligner_Total.png",a1)

(c1 <- ggplot(filter(individual_stats, Enzyme == "PstI-MspI"), aes(Aligner, WGS_intersection))+
    geom_boxplot() +
    labs(x = "Aligner",
         y="SNPs presented in WGS")+
    theme_bw(base_size = 16))
ggsave("PstI_Aligner_WGS.png",c1)

(d1 <- ggplot(filter(individual_stats, Enzyme == "PstI-MspI"), aes(Aligner, ratio))+
    geom_boxplot() +
    labs(x = "Aligner",
         y="SNPs presented in WGS")+
    theme_bw(base_size = 16))
ggsave("PstI_Aligner_Ratio.png",d1)

# Plots without comparison PstI callers
(e1 <- ggplot(filter(individual_stats, Enzyme == "PstI-MspI"), aes(Caller,total_snps))+
    geom_boxplot() +
    labs(x = "Caller",
         y="Total number of SNPs after filtering")+
    theme_bw(base_size = 16)+
    theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust=1)))
ggsave("PstI_Caller_Total.png",e1)

(g1 <- ggplot(filter(individual_stats, Enzyme == "PstI-MspI"), aes(Caller, WGS_intersection))+
    geom_boxplot() +
    labs(x = "Caller",
         y="SNPs presented in WGS")+
    theme_bw(base_size = 16)+
    theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust=1)))
ggsave("PstI_Caller_WGS.png",g1)

(h1 <- ggplot(filter(individual_stats, Enzyme == "PstI-MspI"), aes(Caller, ratio))+
    geom_boxplot() +
    labs(x = "Caller",
         y="SNPs presented in WGS")+
    theme_bw(base_size = 16)+
    theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust=1)))
ggsave("PstI_Caller_Ratio.png",h1)

# Plots without comparison HindIII Aligners
(a2 <- ggplot(filter(individual_stats, Enzyme == "HindIII-NlaIII"), aes(Aligner,total_snps))+
    geom_boxplot() +
    labs(x = "Aligner",
         y="Total number of SNPs after filtering")+
    theme_bw(base_size = 16))
ggsave("HindIII_Aligner_Total.png",a2)

(c2 <- ggplot(filter(individual_stats, Enzyme == "HindIII-NlaIII"), aes(Aligner, WGS_intersection))+
    geom_boxplot() +
    labs(x = "Aligner",
         y="SNPs presented in WGS")+
    theme_bw(base_size = 16))
ggsave("HindIII_Aligner_WGS.png",c2)

(d2 <- ggplot(filter(individual_stats, Enzyme == "HindIII-NlaIII"), aes(Aligner, ratio))+
    geom_boxplot() +
    labs(x = "Aligner",
         y="SNPs presented in WGS")+
    theme_bw(base_size = 16))
ggsave("HindIII_Aligner_Ratio.png",d2)

# Plots without comparison HindIII callers
(e2 <- ggplot(filter(individual_stats, Enzyme == "HindIII-NlaIII"), aes(Caller,total_snps))+
    geom_boxplot() +
    labs(x = "Caller",
         y="Total number of SNPs after filtering")+
    theme_bw(base_size = 16)+
    theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust=1)))
ggsave("HindIII_Caller_Total.png",e2)

(g2 <- ggplot(filter(individual_stats, Enzyme == "HindIII-NlaIII"), aes(Caller, WGS_intersection))+
    geom_boxplot() +
    labs(x = "Aligner",
         y="SNPs presented in WGS")+
    theme_bw(base_size = 16)+
    theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust=1)))
ggsave("HindIII_Caller_WGS.png",g2)

(h2 <- ggplot(filter(individual_stats, Enzyme == "HindIII-NlaIII"), aes(Caller, ratio))+
    geom_boxplot() +
    labs(x = "Caller",
         y="SNPs presented in WGS")+
    theme_bw(base_size = 16)+
    theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust=1)))
ggsave("HindIII_Caller_Ratio.png",h2)

# Aggregative stats for enzyme
(e3 <- ggplot(filter(individual_stats), aes(Enzyme,total_snps))+
    geom_boxplot() +
    labs(x = "Enzyme",
         y="Total number of SNPs after filtering")+
    theme_bw(base_size = 16))
ggsave("All_Total.png",e3)

(g3 <- ggplot(filter(individual_stats), aes(Enzyme, WGS_intersection))+
    geom_boxplot() +
    labs(x = "Enzyme",
         y="SNPs presented in WGS")+
    theme_bw(base_size = 16))
ggsave("All_WGS.png",g3)

(h3 <- ggplot(filter(individual_stats), aes(Enzyme, ratio))+
    geom_boxplot() +
    labs(x = "Enzyme",
         y="SNPs presented in WGS")+
    theme_bw(base_size = 16))
ggsave("All_Ratio.png",h3)


#-------------------------------------------------------------------------

# Pairs plots
full_pairs_test <- mutate(full_pairs, Label_1_new = if_else(Aligner_1 != "Strobealign" & Aligner_2 == "Strobealign" & Enzyme_1 == "HindIII" & Enzyme_2 == "ApeKI",Label_2, Label_1), Label_2_new = if_else(Aligner_1 != "Strobealign" & Aligner_2 == "Strobealign" & Enzyme_1 == "HindIII"& Enzyme_2 == "ApeKI",Label_1, Label_2)) %>% 
  mutate(Label_1_new_new = if_else(Aligner_1 != "Strobealign" & Aligner_2 == "Strobealign" & Enzyme_1 == "PstI" & Enzyme_2 == "HindIII",Label_2_new, Label_1_new), Label_2_new_new = if_else(Aligner_1 != "Strobealign" & Aligner_2 == "Strobealign" & Enzyme_1 == "PstI"& Enzyme_2 == "HindIII",Label_1_new, Label_2_new)) %>% 
  mutate(Label_1_new_new_new = if_else(Aligner_1 != "Strobealign" & Aligner_2 == "Strobealign" & Enzyme_1 == "PstI" & Enzyme_2 == "ApeKI",Label_2_new_new, Label_1_new_new), Label_2_new_new_new = if_else(Aligner_1 != "Strobealign" & Aligner_2 == "Strobealign" & Enzyme_1 == "PstI"& Enzyme_2 == "ApeKI",Label_1_new_new, Label_2_new_new))
ggp <- ggplot(full_pairs_test, aes(Label_1_new_new_new, Label_2_new_new_new)) +  
  geom_tile(aes(fill = R_dosage)) +
  geom_text(aes(label = paste(round(R_dosage, digits = 2), number, sep = "\n"),size = 1))+
  geom_tile(aes(Label_2_new_new_new, Label_1_new_new_new,alpha = NRD))+
  geom_text(aes(Label_2_new_new_new, Label_1_new_new_new,label = round(NRD, digits = 1),size = 1))+
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Label_1",
       y = "Label_2")
  
ggp
ggsave("Pairs_All_Rsquare_NRD.png",ggp, width = 49,height = 44)
ggp <- ggplot(full_pairs, aes(Label_1, Label_2)) +  
  geom_tile(aes(fill = R_dosage)) +
  geom_text(aes(label = round(R_dosage, digits = 2)),size = 3)+
  geom_tile(aes(Label_2, Label_1,alpha = NRD))+
  geom_text(aes(Label_2, Label_1,label = round(NRD, digits = 1),size = 2))+
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggp
ggsave("Pairs_All_Rsquare_NRD.png",ggp, width = 36,height = 32)

ggp1 <- ggplot(filter(full_pairs, Enzyme_1 == "ApeKI" & Enzyme_2 == "ApeKI"), aes(Label_1, Label_2)) +
  geom_tile(aes(fill = R_dosage)) +
  geom_text(aes(label = paste(round(R_dosage, digits = 2), number, sep = "\n"),size = 3))+
  geom_tile(aes(Label_2, Label_1,alpha = NRD))+
  geom_text(aes(Label_2, Label_1,label = paste(round(NRD, digits = 1), TsTv, sep = "\n"),size = 3))+
  theme(text = element_text(size = 18), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_vline(xintercept = 7.5,linewidth=3)+
  geom_vline(xintercept = 14.5,linewidth=3)+
  geom_vline(xintercept = 21.5,linewidth=3)+
  geom_hline(yintercept = 7.5,linewidth=3)+
  geom_hline(yintercept = 14.5,linewidth=3)+
  geom_hline(yintercept = 21.5,linewidth=3)
ggsave("Pairs_ApeKI_Rsquare_NRD.png",ggp1, width = 22,height = 20)

ggp2 <- ggplot(filter(full_pairs, Enzyme_1 == "HindIII" & Enzyme_2 == "HindIII"), aes(Label_1, Label_2)) +
  geom_tile(aes(fill = R_dosage)) +
  geom_text(aes(label = paste(round(R_dosage, digits = 2), number, sep = "\n"),size = 3))+
  geom_tile(aes(Label_2, Label_1,alpha = NRD))+
  geom_text(aes(Label_2, Label_1,label = paste(round(NRD, digits = 1), TsTv, sep = "\n"),size = 3))+
  theme(text = element_text(size = 18), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_vline(xintercept = 7.5,linewidth=3)+
  geom_vline(xintercept = 14.5,linewidth=3)+
  geom_vline(xintercept = 21.5,linewidth=3)+
  geom_hline(yintercept = 7.5,linewidth=3)+
  geom_hline(yintercept = 14.5,linewidth=3)+
  geom_hline(yintercept = 21.5,linewidth=3)
ggsave("Pairs_HindIII_Rsquare_NRD.png",ggp2, width = 22,height = 20)

ggp3 <- ggplot(filter(full_pairs, Enzyme_1 == "PstI" & Enzyme_2 == "PstI"), aes(Label_1, Label_2)) +
  geom_tile(aes(fill = R_dosage)) +
  geom_text(aes(label = paste(round(R_dosage, digits = 2), number, sep = "\n"),size = 3))+
  geom_tile(aes(Label_2, Label_1,alpha = NRD))+
  geom_text(aes(Label_2, Label_1,label = paste(round(NRD, digits = 1), TsTv, sep = "\n"),size = 3))+
  theme(text = element_text(size = 18), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_vline(xintercept = 7.5,linewidth=3)+
  geom_vline(xintercept = 14.5,linewidth=3)+
  geom_vline(xintercept = 21.5,linewidth=3)+
  geom_hline(yintercept = 7.5,linewidth=3)+
  geom_hline(yintercept = 14.5,linewidth=3)+
  geom_hline(yintercept = 21.5,linewidth=3)
ggsave("Pairs_PstI_Rsquare_NRD.png",ggp3, width = 22,height = 20)

full_pairs %>% group_by(Enzyme_pair) %>% summarise(median_intersection=median(number),median_TsTv=median(TsTv))

# SNPs in genes

(i <- ggplot(individual_stats,aes(Enzyme,ratio_genes))+
  geom_boxplot() +
  labs(x = "Enzyme",
       y="Fraction of SNPs in genes")+
  theme_bw(base_size = 16))
ggsave("Genes_ratio_per_enzyme.png",i)

(i1 <- ggplot(filter(individual_stats,total_snps<300000),aes(Enzyme,in_genes))+
    geom_boxplot() +
    labs(x = "Enzyme",
         y="Number of SNPs in genes")+
    theme_bw(base_size = 16))
ggsave("SNPs_in_genes.png",i1)

(i2 <- ggplot(individual_stats,aes(Enzyme,ratio_repeats))+
    geom_boxplot() +
    labs(x = "Enzyme",
         y="Fraction of SNPs in repeats")+
    theme_bw(base_size = 16))
ggsave("Repeats_ratio_per_enzyme.png",i2)

(i3 <- ggplot(filter(individual_stats,total_snps<300000),aes(Enzyme,in_repeats))+
    geom_boxplot() +
    labs(x = "Enzyme",
         y="Number of SNPs in repeats")+
    theme_bw(base_size = 16))
ggsave("SNPs_in_repeats.png",i3)
#-------------------------------------------------------------------------
# SNPs mean r square between different callers
full_pairs_per_enzyme_copy <- full_pairs_per_enzyme
full_pairs_per_enzyme_copy$Caller_1 <- full_pairs_per_enzyme$Caller_2
full_pairs_per_enzyme_copy$Caller_2 <- full_pairs_per_enzyme$Caller_1

full_pairs_per_enzyme_double_caller <- bind_rows(full_pairs_per_enzyme,full_pairs_per_enzyme_copy)
(j1 <- ggplot(full_pairs_per_enzyme_double_caller, aes(Caller_1,R_dosage))+
  geom_boxplot() +
  labs(y = "R_dosage",
       x="Caller")+
  theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, hjust=1)))
ggsave("R_dosage_per_caller.png",j1)

model <- aov(SNPs_ratio ~ Caller,individual_stats)
summary(model)
test <- SNK.test(model, trt = "Caller")
grouping <- select(test$groups, -SNPs_ratio) 
grouping$Caller <- rownames(grouping)
individual_stats_snp <- left_join(individual_stats,grouping, by="Caller")
(e2m <- ggplot(individual_stats_snp, aes(Caller,SNPs_ratio, label=groups))+
    geom_boxplot() +
    geom_text(y= 1.03,alpha=0.2,size=5)+
    labs(x = "Caller",
         y="Amount of SNPs after filtering")+
    theme_bw(base_size = 16)+
    theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust=1)))
ggsave("Caller_Total.png",e2m)


res.aov <- full_pairs_per_enzyme_double_caller %>% kruskal_test(R_dosage ~ Caller_1)
pwc <- full_pairs_per_enzyme_double_caller %>% wilcox_test(R_dosage ~ Caller_1) %>% add_xy_position(x = "Caller_1") %>% arrange(p.adj) %>% head()
j1p <- ggboxplot(full_pairs_per_enzyme_double_caller, x = "Caller_1", y = "R_dosage",xlab="Caller",ylab = "Pearson's r^2 on dosage") +
  stat_pvalue_manual(pwc, hide.ns = TRUE, step.increase = 0.05) +
  labs(
    subtitle = get_test_label(res.aov, detailed = T),
    caption = get_pwc_label(pwc)
  ) + 
  theme_bw(base_size = 16)
ggsave("R_dosage_per_caller_p.png",j1p)

(j2 <- ggplot(full_pairs_per_enzyme_double_caller, aes(Caller_1,NRD))+
  geom_boxplot() +
    labs(y = "NRD",
         x="Caller")+
  theme_bw(base_size = 16))
ggsave("NRD_dosage_per_caller.png",j2)
(j3 <- ggplot(full_pairs_per_enzyme_double_caller, aes(Caller_1,number))+
  geom_boxplot() +
    labs(y = "SNPs",
         x="Caller")+
  theme_bw(base_size = 16))
ggsave("SNPs_per_caller.png",j3)
(j4 <- ggplot(full_pairs_per_enzyme_double_caller, aes(Caller_1,TsTv))+
  geom_boxplot() +
    labs(y = "TsTv",
         x="Caller")+
  theme_bw(base_size = 16))
ggsave("TsTv_per_caller.png",j4)
full_pairs_per_enzyme_double_caller %>% group_by(Caller_1) %>% tally()

#-------------------------------------------------------------------------
# SNPs mean r square between different aligners
full_pairs_per_enzyme_copy <- full_pairs_per_enzyme
full_pairs_per_enzyme_copy$Aligner_1 <- full_pairs_per_enzyme$Aligner_2
full_pairs_per_enzyme_copy$Aligner_2 <- full_pairs_per_enzyme$Aligner_1

full_pairs_per_enzyme_double_align <- bind_rows(full_pairs_per_enzyme,full_pairs_per_enzyme_copy)
(k1 <- ggplot(full_pairs_per_enzyme_double_align, aes(Aligner_1,R_dosage))+
    geom_boxplot() +
    labs(y = "R_dosage",
         x="Aligner")+
    theme_bw(base_size = 16))
ggsave("R_dosage_per_aligner.png",k1)

res.aov <- full_pairs_per_enzyme_double_align %>% kruskal_test(R_dosage ~ Aligner_1)
pwc <- full_pairs_per_enzyme_double_align %>% wilcox_test(R_dosage ~ Aligner_1) %>% add_xy_position(x = "Aligner_1")
k1p <- ggboxplot(full_pairs_per_enzyme_double_align, x = "Aligner_1", y = "R_dosage",xlab="Aligner",ylab = "Pearson's r^2 on dosage") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = T),
    caption = get_pwc_label(pwc)
  ) + 
  theme_bw(base_size = 16)
ggsave("R_dosage_per_aligner_p.png",k1p)

model <- aov(R_dosage ~ Aligner_1,full_pairs_per_enzyme_double_align)

(k2 <- ggplot(full_pairs_per_enzyme_double_align, aes(Aligner_1,NRD))+
    geom_boxplot() +
    labs(y = "NRD",
         x="Aligner")+
    theme_bw(base_size = 16))
ggsave("NRD_dosage_per_aligner.png",k2)
(k3 <- ggplot(full_pairs_per_enzyme_double_align, aes(Aligner_1,number))+
    geom_boxplot() +
    labs(y = "SNPs",
         x="Aligner")+
    theme_bw(base_size = 16))
ggsave("SNPs_per_aligner.png",k3)
(k4 <- ggplot(full_pairs_per_enzyme_double_align, aes(Aligner_1,TsTv))+
    geom_boxplot() +
    labs(y = "TsTv",
         x="Aligner")+
    theme_bw(base_size = 16))
ggsave("TsTv_per_aligner.png",k4)
full_pairs_per_enzyme_double_align %>% group_by(Aligner_1) %>% tally()
#-------------------------------------------------------------------------
# SNPs mean r square between different enzyme
full_pairs_per_enzyme_copy <- full_pairs_per_enzyme
#full_pairs_per_enzyme_copy$Enzyme_1 <- full_pairs_per_enzyme$Enzyme_2
#full_pairs_per_enzyme_copy$Enzyme_2 <- full_pairs_per_enzyme$Enzyme_1

#full_pairs_per_enzyme_double_enzyme <- bind_rows(full_pairs_per_enzyme,full_pairs_per_enzyme_copy)
full_pairs_per_enzyme$Enzyme_1 <- factor(full_pairs_per_enzyme$Enzyme_1, labels = c("ApeKI","HindIII-NlaIII","PstI-MspI"))
(k1 <- ggplot(full_pairs_per_enzyme, aes(Enzyme_1,R_dosage))+
    geom_boxplot() +
    labs(y = "R_dosage",
         x="Enzyme set")+
    theme_bw(base_size = 16))
ggsave("R_dosage_per_enzyme.png",k1)

res.aov <- full_pairs_per_enzyme %>% kruskal_test(R_dosage ~ Enzyme_1)
pwc <- full_pairs_per_enzyme %>% wilcox_test(R_dosage ~ Enzyme_1) %>% add_xy_position(x = "Enzyme_1")
k1p <- ggboxplot(full_pairs_per_enzyme, x = "Enzyme_1", y = "R_dosage",xlab="Enzyme",ylab = "Pearson's r^2 on dosage") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = T),
    caption = get_pwc_label(pwc)
  ) + 
  theme_bw(base_size = 16)
ggsave("R_dosage_per_enzyme_p.png",k1p)

(k2 <- ggplot(full_pairs_per_enzyme, aes(Enzyme_1,NRD))+
    geom_boxplot() +
    labs(y = "NRD",
         x="Enzymes set")+
    theme_bw(base_size = 16))
ggsave("NRD_dosage_per_enzyme.png",k2)
(k3 <- ggplot(full_pairs_per_enzyme, aes(Enzyme_1,number))+
    geom_boxplot() +
    labs(y = "SNPs",
         x="Aligner")+
    theme_bw(base_size = 16))
ggsave("SNPs_per_enzyme.png",k3)
(k4 <- ggplot(full_pairs_per_enzyme, aes(Enzyme_1,TsTv))+
    geom_boxplot() +
    labs(y = "TsTv",
         x="Aligner")+
    theme_bw(base_size = 16))
ggsave("TsTv_per_enzyme.png",k4)
full_pairs_per_enzyme_double_enzyme %>% group_by(Enzyme_1) %>% tally()
#-------------------------------------------------------------------------
df <- data_frame(SNPs = c(119007, 35531, 21648, 20930, 24117, 19110, 853), callers = c(1,2,3,4,5,6,7))
df$callers <- as.factor(df$callers)
snps_per_callers <- ggplot(df, aes(x=callers,y=SNPs))+
  geom_col()+
  theme_bw(base_size = 16)
ggsave("SNPs_per_caller_from_venn.png",snps_per_callers)
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
TP_TN_ApeKI_BWA <- read.csv("~/Yandex.Disk.localized/Skoltech/DivSoy/Benchmarking/Results/Benchmarking_18M/Depth_5/TP_TN/TP_TN_data_ApeKI_BWA.csv")
TP_TN_ApeKI_BWA <-  mutate(TP_TN_ApeKI_BWA, Enzyme_set = "ApeKI", Aligner = "BWA")
TP_TN_HindIII_BWA <- read.csv("~/Yandex.Disk.localized/Skoltech/DivSoy/Benchmarking/Results/Benchmarking_18M/Depth_5/TP_TN/TP_TN_data_HindIII_BWA.csv")
TP_TN_HindIII_BWA <-  mutate(TP_TN_HindIII_BWA, Enzyme_set = "HindIII-NlaIII", Aligner = "BWA")
TP_TN_PstI_BWA <- read.csv("~/Yandex.Disk.localized/Skoltech/DivSoy/Benchmarking/Results/Benchmarking_18M/Depth_5/TP_TN/TP_TN_data_PstI_BWA.csv")
TP_TN_PstI_BWA <-  mutate(TP_TN_PstI_BWA, Enzyme_set = "PstI-MspI", Aligner = "BWA")
TP_TN_ApeKI_BBMap <- read.csv("~/Yandex.Disk.localized/Skoltech/DivSoy/Benchmarking/Results/Benchmarking_18M/Depth_5/TP_TN/TP_TN_data_ApeKI_BBMap.csv")
TP_TN_ApeKI_BBMap <-  mutate(TP_TN_ApeKI_BBMap, Enzyme_set = "ApeKI", Aligner = "BBMap")
TP_TN_HindIII_BBMap <- read.csv("~/Yandex.Disk.localized/Skoltech/DivSoy/Benchmarking/Results/Benchmarking_18M/Depth_5/TP_TN/TP_TN_data_HindIII_BBMap.csv")
TP_TN_HindIII_BBMap <-  mutate(TP_TN_HindIII_BBMap, Enzyme_set = "HindIII-NlaIII", Aligner = "BBMap")
TP_TN_PstI_BBMap <- read.csv("~/Yandex.Disk.localized/Skoltech/DivSoy/Benchmarking/Results/Benchmarking_18M/Depth_5/TP_TN/TP_TN_data_PstI_BBMap.csv")
TP_TN_PstI_BBMap <-  mutate(TP_TN_PstI_BBMap, Enzyme_set = "PstI-MspI", Aligner = "BBMap")
TP_TN_ApeKI_Bowtie2 <- read.csv("~/Yandex.Disk.localized/Skoltech/DivSoy/Benchmarking/Results/Benchmarking_18M/Depth_5/TP_TN/TP_TN_data_ApeKI_Bowtie2.csv")
TP_TN_ApeKI_Bowtie2 <-  mutate(TP_TN_ApeKI_Bowtie2, Enzyme_set = "ApeKI", Aligner = "Bowtie2")
TP_TN_HindIII_Bowtie2 <- read.csv("~/Yandex.Disk.localized/Skoltech/DivSoy/Benchmarking/Results/Benchmarking_18M/Depth_5/TP_TN/TP_TN_data_HindIII_Bowtie2.csv")
TP_TN_HindIII_Bowtie2 <-  mutate(TP_TN_HindIII_Bowtie2, Enzyme_set = "HindIII-NlaIII", Aligner = "Bowtie2")
TP_TN_PstI_Bowtie2 <- read.csv("~/Yandex.Disk.localized/Skoltech/DivSoy/Benchmarking/Results/Benchmarking_18M/Depth_5/TP_TN/TP_TN_data_PstI_Bowtie2.csv")
TP_TN_PstI_Bowtie2 <-  mutate(TP_TN_PstI_Bowtie2, Enzyme_set = "PstI-MspI", Aligner = "Bowtie2")
TP_TN_ApeKI_str <- read.csv("~/Yandex.Disk.localized/Skoltech/DivSoy/Benchmarking/Results/Benchmarking_18M/Depth_5/TP_TN/TP_TN_data_ApeKI_str.csv")
TP_TN_ApeKI_str <-  mutate(TP_TN_ApeKI_str, Enzyme_set = "ApeKI", Aligner = "Strobealign")
TP_TN_HindIII_str <- read.csv("~/Yandex.Disk.localized/Skoltech/DivSoy/Benchmarking/Results/Benchmarking_18M/Depth_5/TP_TN/TP_TN_data_HindIII_str.csv")
TP_TN_HindIII_str <-  mutate(TP_TN_HindIII_str, Enzyme_set = "HindIII-NlaIII", Aligner = "Strobealign")
TP_TN_PstI_str <- read.csv("~/Yandex.Disk.localized/Skoltech/DivSoy/Benchmarking/Results/Benchmarking_18M/Depth_5/TP_TN/TP_TN_data_PstI_str.csv")
TP_TN_PstI_str <-  mutate(TP_TN_PstI_str, Enzyme_set = "PstI-MspI", Aligner = "Strobealign")
TP_TN_data <- bind_rows(TP_TN_ApeKI_BWA,TP_TN_HindIII_BWA,TP_TN_PstI_BWA,TP_TN_ApeKI_BBMap,TP_TN_HindIII_BBMap,TP_TN_PstI_BBMap,TP_TN_ApeKI_Bowtie2,TP_TN_HindIII_Bowtie2,TP_TN_PstI_Bowtie2,TP_TN_ApeKI_str,TP_TN_HindIII_str,TP_TN_PstI_str)

TP_TN_data <- TP_TN_data %>%
  mutate(
    F1 = 2 * Precision * Recall / (Precision + Recall)
    # , F1_alt = 2 * TP / (2*TP + FP + FN)   # alternate formula
  )
colnames(TP_TN_data)[1] <- "Caller"
View(TP_TN_data)
mean(filter(TP_TN_data,Caller == "DeepVariant")$Precision)
mean(filter(TP_TN_data,Caller == "DeepVariant")$Recall)
mean(filter(TP_TN_data,Caller == "DeepVariant")$FPR)
mean(filter(TP_TN_data,Caller == "freebayes")$Precision)
mean(filter(TP_TN_data,Caller == "freebayes")$FPR)
write.csv(TP_TN_data, "TP_TN_data.csv")
(precision1 <- ggplot(TP_TN_data, aes(x=Caller,y=Precision,fill=Enzyme_set))+
    geom_boxplot()+
    theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ylim(0,1))
ggsave("precision_enzyme_set.png",precision1)
(precision2 <- ggplot(TP_TN_data, aes(x=Caller,y=Precision,fill=Aligner))+
    geom_boxplot()+
    theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ylim(0,1))
ggsave("precision_aligner.png",precision2)

(recall1 <- ggplot(TP_TN_data, aes(x=Caller,y=Recall,fill=Enzyme_set))+
    geom_boxplot()+
    theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ylim(0,1))
ggsave("recall_enzyme_set.png",recall1)
(recall2 <- ggplot(TP_TN_data, aes(x=Caller,y=Recall,fill=Aligner))+
    geom_boxplot()+
    theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ylim(0,1))
ggsave("recall_aligner.png",recall2)

(fpr1 <- ggplot(TP_TN_data, aes(x=Caller,y=FPR,fill=Enzyme_set))+
    geom_boxplot()+
    theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ylim(0,1))
ggsave("FPR_enzyme_set.png",fpr1)
(fpr2 <- ggplot(TP_TN_data, aes(x=Caller,y=FPR,fill=Aligner))+
    geom_boxplot()+
    theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ylim(0,1))
ggsave("FPR_aligner.png",fpr2)


(f1_1 <- ggplot(TP_TN_data, aes(x=Caller,y=F1,fill=Enzyme_set))+
    geom_boxplot()+
    theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ylim(0,1))
ggsave("F1_enzyme_set.png",fpr1)
(f1_2 <- ggplot(TP_TN_data, aes(x=Caller,y=F1,fill=Aligner))+
    geom_boxplot()+
    theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ylim(0,1))
ggsave("F1_aligner.png",fpr2)

#-------------------------------------------------------------------------
# ANOVA tests
model_ratio_enzyme <- aov(ratio ~ Enzyme,individual_stats)
summary(model_ratio_enzyme)
shapiro.test(residuals(model_aov)) 

# Проверка гомогенности дисперсии
bartlett.test( residuals(model_aov) ~ individual_stats$Enzyme) # https://ru.wikipedia.org/wiki/Критерий_Бартлетта
leveneTest( residuals(model_aov) ~ individual_stats$Enzyme) # https://en.wikipedia.org/wiki/Levene%27s_test

model_total_snps_enzyme <- aov(total_snps ~ Enzyme,individual_stats)
summary(model_total_snps_enzyme)
model_ratio_genes_enzyme <- aov(ratio_genes ~ Enzyme,individual_stats)
summary(model_ratio_genes_enzyme)
model_in_genes_enzyme <- aov(in_genes ~ Enzyme,individual_stats)
summary(model_in_genes_enzyme)
t.test(filter(individual_stats, Enzyme == "ApeKI")$total_snps,filter(individual_stats, Enzyme == "HindIII-NlaIII")$total_snps)

individual_stats$total_snps <- as.double(individual_stats$total_snps)
my_comparisons <- list( c("ApeKI", "HindIII-NlaIII"), c("ApeKI", "PstI-MspI"), c("HindIII-NlaIII", "PstI-MspI") )
e3p <- ggboxplot(individual_stats, x = "Enzyme", y = "total_snps",ylab = "Total number of SNPs after filtering") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + 
  theme_bw(base_size = 16) 
ggsave("All_Total_p.png",e3p)
ggexport(e3p,filename="All_Total_p.png",res=30)  





i1p <- ggboxplot(filter(individual_stats,total_snps<300000), x = "Enzyme", y = "in_genes",ylab = "Number of SNPs in genes") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + 
  theme_bw(base_size = 16) 
ggsave("SNPs_in_genes_p.png",i1p)

(i1 <- ggplot(filter(individual_stats,total_snps<300000),aes(Enzyme,in_genes))+
    geom_boxplot() +
    labs(x = "Enzyme",
         y="Number of SNPs in genes")+
    theme_bw(base_size = 16))
ggsave("SNPs_in_genes.png",i1)

res.aov <- individual_stats %>% anova_test(total_snps ~ Enzyme)
pwc <- individual_stats %>% tukey_hsd(total_snps ~ Enzyme) %>% add_xy_position(x = "Enzyme")
e3p <- ggboxplot(individual_stats, x = "Enzyme", y = "total_snps",ylab = "Total number of SNPs after filtering") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = T),
    caption = get_pwc_label(pwc)
  ) + 
  theme_bw(base_size = 16)
ggsave("All_Total_p.png",e3p)

res.aov <- individual_stats %>% mutate(ratio_trans=asin(sqrt(ratio))) %>% anova_test(ratio_trans ~ Enzyme)

pwc <- individual_stats %>% mutate(ratio_trans=asin(sqrt(ratio))) %>% tukey_hsd(ratio_trans ~ Enzyme) %>% add_xy_position(x = "Enzyme") %>% mutate(y.position = y.position -0.25)
h3p <- ggboxplot(individual_stats, x = "Enzyme", y = "ratio",ylab = "SNPs presented in WGS") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = T),
    caption = get_pwc_label(pwc)
  ) + 
  theme_bw(base_size = 16)
ggsave("All_Ratio_p_transformed.png",h3p)

 i1p <- ggboxplot(filter(individual_stats,total_snps<300000), x = "Enzyme", y = "in_genes",ylab = "Number of SNPs in genes") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + 
  theme_bw(base_size = 16) 
ggsave("SNPs_in_genes_p.png",i1p)

res.aov <- individual_stats %>% anova_test(in_genes ~ Enzyme)
pwc <- individual_stats %>% tukey_hsd(in_genes ~ Enzyme) %>% add_xy_position(x = "Enzyme")
i1p <- ggboxplot(filter(individual_stats,total_snps<300000), x = "Enzyme", y = "in_genes",ylab = "Number of SNPs in genes") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = T),
    caption = get_pwc_label(pwc)
  ) + 
  theme_bw(base_size = 16)
ggsave("SNPs_in_genes_p.png",i1p)


res.aov <- individual_stats %>% mutate(ratio_genes_trans=asin(sqrt(ratio_genes)))%>% anova_test(ratio_genes_trans ~ Enzyme)
pwc <- individual_stats %>% mutate(ratio_genes_trans=asin(sqrt(ratio_genes)))%>% tukey_hsd(ratio_genes_trans ~ Enzyme) %>% add_xy_position(x = "Enzyme") %>% mutate(y.position = y.position -0.3)
ip <- ggboxplot(individual_stats, x = "Enzyme", y = "ratio_genes",ylab = "Fraction of SNPs in genes") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = T),
    caption = get_pwc_label(pwc)
  ) + 
  theme_bw(base_size = 16) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1))
ggsave("Genes_ratio_per_enzyme_p_transformed.png",ip)

res.aov <- individual_stats %>% anova_test(ratio ~ Enzyme+Aligner)
pwc <- individual_stats %>% tukey_hsd(ratio ~ Aligner + Enzyme) %>% add_xy_position(x = "Enzyme")
ip <- ggboxplot(individual_stats, x = "Aligner", y = "ratio",facet.by = "Enzyme",ylab = "Fraction of SNPs in genes") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = T),
    caption = get_pwc_label(pwc)
  ) + 
  theme_bw(base_size = 16)
ggsave("Genes_ratio_per_enzyme_p.png",ip)





res.aov <- individual_stats %>% mutate(ratio_trans=asin(sqrt(ratio))) %>% anova_test(ratio_trans ~ Aligner)
pwc <- individual_stats %>% mutate(ratio_trans=asin(sqrt(ratio))) %>% tukey_hsd(ratio_trans ~ Aligner) %>% add_xy_position(x = "Aligner")
d5p <- ggboxplot(individual_stats, x = "Aligner", y = "ratio",ylab = "SNPs presented in WGS") +
  stat_pvalue_manual(pwc,hide.ns = T) +
  labs(
    subtitle = get_test_label(res.aov, detailed = T),
    caption = get_pwc_label(pwc)
  ) + 
  theme_bw(base_size = 16)
ggsave("Aligner_Ratio_transformed.png",d5p)

model <- lm(total_snps ~ Enzyme+Aligner, data = individual_stats)
res.aov <- individual_stats %>% anova_test(total_snps ~ Enzyme+Aligner)
individual_stats %>%
  group_by(Enzyme) %>%
  anova_test(total_snps ~ Aligner, error = model)
pwc <- individual_stats %>% 
  emmeans_test(
    total_snps ~ Aligner, p.adjust.method = "bonferroni",
    model = model
  )
res.aov <- individual_stats %>% mutate(SNPs_ratio_trans=asin(sqrt(SNPs_ratio))) %>% anova_test(SNPs_ratio_trans ~ Aligner)
pwc <- individual_stats %>% mutate(SNPs_ratio_trans=asin(sqrt(SNPs_ratio))) %>% tukey_hsd(SNPs_ratio_trans ~ Aligner) %>% add_xy_position(x = "Aligner")
a5p <- ggboxplot(individual_stats, x = "Aligner", y = "SNPs_ratio", ylab = "Amount of SNPs after filtering") +
  stat_pvalue_manual(pwc,hide.ns = T) +
  labs(
    subtitle = get_test_label(res.aov, detailed = T),
    caption = get_pwc_label(pwc)
  ) + 
  theme_bw(base_size = 16)
ggsave("Aligner_Total_transformed.png",a5p)

(a5 <- ggplot(filter(individual_stats), aes(Enzyme,total_snps,fill=Aligner))+
    geom_boxplot() +
    labs(x = "Enzyme set",
         y="Total number of SNPs")+
    theme_bw(base_size = 16))

# Precision, Recall, FPR plots
TP_TN_data <- mutate(TP_TN_data,Precision_trans = asin(sqrt(Precision)))
test <- HSD.test(aov(Precision_trans ~ Caller,TP_TN_data), trt = "Caller")
grouping <- select(test$groups, -Precision_trans) 
grouping$Caller <- rownames(grouping)
TP_TN_data_prec <- left_join(TP_TN_data,grouping, by="Caller")
(precision1 <- ggplot(TP_TN_data_prec, aes(x=Caller,y=Precision, label=groups))+
    geom_boxplot()+
    geom_text(y= 1.03,alpha=0.2,size=5)+
    theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ylim(0,1.05))
ggsave("precision_transformed.png",precision1)

TP_TN_data <- mutate(TP_TN_data,Recall_trans = asin(sqrt(Recall)))
test <- HSD.test(aov(Recall_trans ~ Caller,TP_TN_data), trt = "Caller")
grouping <- select(test$groups, -Recall_trans) 
grouping$Caller <- rownames(grouping)
TP_TN_data_recall <- left_join(TP_TN_data,grouping, by="Caller")
(recall1 <- ggplot(TP_TN_data_recall, aes(x=Caller,y=Recall, label=groups))+
    geom_boxplot()+
    geom_text(y= 1.03,alpha=0.2,size=5)+
    theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ylim(0,1.05))
ggsave("recall_transformed.png",recall1)

summary(aov(Precision ~ Caller,TP_TN_data))

TP_TN_data <- mutate(TP_TN_data,FPR_trans = asin(sqrt(FPR)))
test <- HSD.test(aov(FPR_trans ~ Caller,TP_TN_data), trt = "Caller")
grouping <- select(test$groups, -FPR_trans) 
grouping$Caller <- rownames(grouping)
TP_TN_data_FPR <- left_join(TP_TN_data,grouping, by="Caller")
(FPR1 <- ggplot(TP_TN_data_FPR, aes(x=Caller,y=FPR, label=groups))+
    geom_boxplot()+
    geom_text(y= 1.03,alpha=0.2,size=5)+
    theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ylim(0,1.05)+
    labs(y="FDR"))
ggsave("FPR_tranformed.png",FPR1)

TP_TN_data <- mutate(TP_TN_data,F1_trans = asin(sqrt(F1)))
test <- HSD.test(aov(F1_trans ~ Caller,TP_TN_data), trt = "Caller")
grouping <- select(test$groups, -F1_trans) 
grouping$Caller <- rownames(grouping)
TP_TN_data_F1 <- left_join(TP_TN_data,grouping, by="Caller")
(F1_trans <- ggplot(TP_TN_data_F1, aes(x=Caller,y=F1, label=groups))+
    geom_boxplot()+
    geom_text(y= 1.03,alpha=0.2,size=5)+
    theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ylim(0,1.05)+
    labs(y="F1"))
ggsave("F1_tranformed.png",F1_trans)

res.aov <- TP_TN_data %>% anova_test(Precision ~ Caller)
pwc <- TP_TN_data %>% tukey_hsd(Precision ~ Caller) %>% add_xy_position(x = "Caller")
SNK
precision1_p <- ggboxplot(TP_TN_data, x = "Caller", y = "Precision") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = T),
    caption = get_pwc_label(pwc)
  ) + 
  theme_bw(base_size = 16)
ggsave("All_Total_p.png",e3p)

(recall1 <- ggplot(TP_TN_data, aes(x=Caller,y=Recall))+
    geom_boxplot()+
    theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ylim(0,1))
ggsave("recall_enzyme_set.png",recall1)

(fpr <- ggplot(TP_TN_data, aes(x=Caller,y=FPR))+
    geom_boxplot()+
    theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ylim(0,1))
ggsave("FPR_enzyme_set.png",fpr1)

individual_stats <- individual_stats %>% group_by(Enzyme) %>% mutate(SNPs_ratio = total_snps/max(total_snps)) %>% ungroup()
individual_stats %>% group_by(Caller) %>% summarise(mean(ratio))
individual_stats %>% group_by(Caller) %>% summarise(mean(SNPs_ratio))

individual_stats <- mutate(individual_stats, SNPs_ratio_trans = asin(sqrt(SNPs_ratio)))

model <- aov(SNPs_ratio_trans ~ Caller,individual_stats)
summary(model)

test <- HSD.test(model, trt = "Caller")
grouping <- select(test$groups, -SNPs_ratio_trans) 
grouping$Caller <- rownames(grouping)
individual_stats_snp <- left_join(individual_stats,grouping, by="Caller")
(e2m <- ggplot(individual_stats_snp, aes(Caller,SNPs_ratio, label=groups))+
    geom_boxplot() +
    geom_text(y= 1.03,alpha=0.2,size=5)+
    labs(x = "Caller",
         y="Amount of SNPs after filtering")+
    theme_bw(base_size = 16)+
    theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust=1)))
ggsave("Caller_Total_transformed.png",e2m)

individual_stats <- mutate(individual_stats, ratio_trans = asin(sqrt(ratio)))
model <- aov(ratio_trans ~ Caller,individual_stats)
summary(model)
test <- HSD.test(model, trt = "Caller")
grouping <- select(test$groups, -ratio_trans) 
grouping$Caller <- rownames(grouping)
individual_stats_ratio <- left_join(individual_stats,grouping, by="Caller")
(e2r <- ggplot(individual_stats_ratio, aes(Caller,ratio, label=groups))+
    geom_boxplot() +
    geom_text(y= 0.82,alpha=0.2,size=5)+
    labs(x = "Caller",
         y="SNPs presented in WGS")+
    theme_bw(base_size = 16)+
    theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust=1)))
ggsave("Caller_Ratio_transformed.png",e2r)



(e2m <- ggplot(filter(individual_stats), aes(Caller,SNPs_ratio, label=groups))+
    geom_boxplot() +
    labs(x = "Caller",
         y="Total number of SNPs after filtering")+
    theme_bw(base_size = 16)+
    theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust=1)))
ggsave("HindIII_Caller_Total.png",e2)


