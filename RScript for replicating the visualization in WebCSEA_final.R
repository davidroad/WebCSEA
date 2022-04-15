##Rscript for replicating the visualization in WebCSEA (https://bioinfo.uth.edu/webcsea/).
##Here we use the Example gene list in the WebCSEA from Hepatocyte signature genes
library(ggplot2)
library(ggrepel)
library(tidyverse)
input_path <- "/PATH/TO/WebCSEA/DOWNLOAD/RESULTS/"
output_path <- "/PATH/TO/WORK/DIRECTORY/"
run_name <- "NAME_YOUR_PROJECT"

setwd(output_path)
##please find the WebCSEA.all_pvalue.txt file from the Download results zip folder.
res <- read.delim(paste0(input_path,"WebCSEA_All_tissue_cell_type.txt"),as.is=T,check.names=F)
label1 = res %>% arrange(desc(input_list_combined_log10p)) %>% head(5) %>%
    mutate(delabel = Tissue_cell_type_name)%>% select(Tissue_cell_type_name, delabel)
res = res %>% left_join(label1)

p1_clean <- ggplot(res, aes(x=Organ_system, y=input_list_combined_log10p,col=Organ_system))  + ggtitle("Organ system") + 
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res$input_list_combined_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_clean_Organ_system.pdf"),12,8)
print(p1_clean)
dev.off()

p1 <- ggplot(res, aes(x=Organ_system, y=input_list_combined_log10p,col=Organ_system,label = delabel))+ geom_text_repel(size = 3) + ggtitle("Organ system") + 
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res$input_list_combined_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_Organ_system.pdf"),12,8)
print(p1)
dev.off()

##top20 tissue cell_types result##
res_ordered <- res[order(res$input_list_combined_log10p,decreasing =T),]
top_20_tissue <- unique(res_ordered$Tissue)[1:20]
tissue_idx <- {}
for (i in 1:length(top_20_tissue)){
	tissue_idx <- c(tissue_idx, which(res_ordered$Tissue==top_20_tissue[i]))
}
res_ordered_top20 <- res_ordered[tissue_idx,]
res_ordered_top20$Tissue <- ordered(res_ordered_top20$Tissue,levels = unique(res_ordered_top20[order(res_ordered_top20$input_list_combined_log10p,decreasing=T),4]))

label2 = res_ordered_top20  %>% select(-delabel)%>% arrange(desc(input_list_combined_log10p)) %>% head(5) %>%
    mutate(delabel = Specific_cell_type)%>% select(Tissue_cell_type_name, delabel)
res_ordered_top20 = res_ordered_top20  %>% select(-delabel)%>% left_join(label2)

p2_clean <- ggplot(res_ordered_top20, aes(x=Tissue, y=input_list_combined_log10p,col=Tissue)) + ggtitle("Top20 tissues") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res_ordered$input_list_combined_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_clean_Top20_tissues.pdf"),12,8)
print(p2_clean)
dev.off()

p2 <- ggplot(res_ordered_top20, aes(x=Tissue, y=input_list_combined_log10p,col=Tissue, label = delabel)) +geom_text_repel(size = 3) + ggtitle("Top20 tissues") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res_ordered$input_list_combined_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_Top20_tissues.pdf"),12,8)
print(p2)
dev.off()
##Adult and Fetal res
res_Adult <- res[res$Developmental_stage=="Adult",]
res_Fetal <- res[res$Developmental_stage=="Fetal",]

label3 = res_Adult %>% select(-delabel) %>% arrange(desc(input_list_combined_log10p)) %>% head(5) %>%
    mutate(delabel = Tissue_cell_type_name)%>% select(Tissue_cell_type_name, delabel)
res_Adult = res_Adult  %>% select(-delabel)%>% left_join(label3)

label4 = res_Fetal %>% select(-delabel) %>% arrange(desc(input_list_combined_log10p)) %>% head(5) %>%
    mutate(delabel = Tissue_cell_type_name)%>% select(Tissue_cell_type_name, delabel)
res_Fetal = res_Fetal  %>% select(-delabel)%>% left_join(label4)

p3_clean <- ggplot(res_Adult, aes(x=Organ_system, y=input_list_combined_log10p,col=Organ_system)) + ggtitle("Adult Organ system") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res_Adult$input_list_combined_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_clean_Adult_Organ_system.pdf"),12,8)
print(p3_clean)
dev.off()

p3 <- ggplot(res_Adult, aes(x=Organ_system, y=input_list_combined_log10p,col=Organ_system, label = delabel)) + geom_text_repel(size = 3)+ ggtitle("Adult Organ system") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res_Adult$input_list_combined_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_Adult_Organ_system.pdf"),12,8)
print(p3)
dev.off()


p4_clean <- ggplot(res_Fetal, aes(x=Organ_system, y=input_list_combined_log10p,col=Organ_system)) + ggtitle("Fetal Organ system") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res_Fetal$input_list_combined_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_clean_Fetal_Organ_system.pdf"),12,8)
print(p4_clean)
dev.off()

p4 <- ggplot(res_Fetal, aes(x=Organ_system, y=input_list_combined_log10p,col=Organ_system, label = delabel))+ geom_text_repel(size = 3) + ggtitle("Fetal Organ system") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res_Fetal$input_list_combined_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_Fetal_Organ_system.pdf"),12,8)
print(p4)
dev.off()
##top 20 general cell types##


res_ordered <- res[order(res$input_list_combined_log10p,decreasing =T),]
top_20_cell_type <- unique(res_ordered$General_cell_type)[1:20]
cell_type_idx <- {}
for (i in 1:length(top_20_cell_type)){
	cell_type_idx <- c(cell_type_idx, which(res_ordered$General_cell_type==top_20_cell_type[i]))
}
res_ordered_top20_CT <- res_ordered[cell_type_idx,]
res_ordered_top20_CT$General_cell_type <- ordered(res_ordered_top20_CT$General_cell_type,levels = unique(res_ordered_top20_CT[order(res_ordered_top20_CT$input_list_combined_log10p,decreasing=T),6]))
label5 = res_ordered_top20_CT  %>% select(-delabel)%>% arrange(desc(input_list_combined_log10p)) %>% head(5) %>%
    mutate(delabel = Specific_cell_type)%>% select(Tissue_cell_type_name, delabel)
res_ordered_top20_CT = res_ordered_top20_CT  %>% select(-delabel)%>% left_join(label5)

p5_clean <- ggplot(res_ordered_top20_CT, aes(x=General_cell_type, y=input_list_combined_log10p,col=General_cell_type)) + ggtitle("Top20 general cell types") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res_ordered_top20_CT$input_list_combined_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_clean_Top20_general_cell_type.pdf"),12,8)
print(p5_clean)
dev.off()

p5 <- ggplot(res_ordered_top20_CT, aes(x=General_cell_type, y=input_list_combined_log10p,col=General_cell_type, label = delabel))+ geom_text_repel(size = 3) + ggtitle("Top20 general cell types") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res_ordered_top20_CT$input_list_combined_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_Top20_general_cell_type.pdf"),12,8)
print(p5)
dev.off()
##
##generate the raw pvalue plots
res <- data.frame(res, input_list_raw_log10p = -log10(unlist(res$input_list_raw_p)), check.names=F)
label6 = res %>% arrange(desc(input_list_raw_log10p)) %>% head(5) %>%
    mutate(delabel = Tissue_cell_type_name)%>% select(Tissue_cell_type_name, delabel)
res = res  %>% left_join(label6)

p6_clean <- ggplot(res, aes(x=Organ_system, y=input_list_raw_log10p,col=Organ_system)) + ggtitle("Organ system") + 
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res$input_list_raw_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_raw_pvalue_clean_Organ_system.pdf"),12,8)
print(p6_clean)
dev.off()

p6 <- ggplot(res, aes(x=Organ_system, y=input_list_raw_log10p,col=Organ_system, label = delabel)) + geom_text_repel(size = 3) + ggtitle("Organ system") + 
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res$input_list_raw_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_raw_pvalue_Organ_system.pdf"),12,8)
print(p6)
dev.off()

##top20 tissue cell_types result##
res_ordered <- res[order(res$input_list_raw_log10p,decreasing =T),]
top_20_tissue <- unique(res_ordered$Tissue)[1:20]
tissue_idx <- {}
for (i in 1:length(top_20_tissue)){
	tissue_idx <- c(tissue_idx, which(res_ordered$Tissue==top_20_tissue[i]))
}
res_ordered_top20 <- res_ordered[tissue_idx,]
res_ordered_top20$Tissue <- ordered(res_ordered_top20$Tissue,levels = unique(res_ordered_top20[order(res_ordered_top20$input_list_raw_log10p,decreasing=T),4]))
label7 = res_ordered_top20  %>% select(-delabel)%>% arrange(desc(input_list_raw_log10p)) %>% head(5) %>%
    mutate(delabel = Specific_cell_type)%>% select(Tissue_cell_type_name, delabel)
res_ordered_top20 = res_ordered_top20  %>% select(-delabel)%>% left_join(label7)

p7_clean <- ggplot(res_ordered_top20, aes(x=Tissue, y=input_list_raw_log10p,col=Tissue)) + ggtitle("Top20 tissues") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res_ordered_top20$input_list_raw_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_raw_pvalue_clean_Top20_tissues.pdf"),12,8)
print(p7_clean)
dev.off()

p7 <- ggplot(res_ordered_top20, aes(x=Tissue, y=input_list_raw_log10p,col=Tissue, label = delabel)) + geom_text_repel(size = 3) + ggtitle("Top20 tissues") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res_ordered_top20$input_list_raw_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_raw_pvalue_Top20_tissues.pdf"),12,8)
print(p7)
dev.off()
##Adult and Fetal res
res_Adult <- res[res$Developmental_stage=="Adult",]
res_Fetal <- res[res$Developmental_stage=="Fetal",]

label8 = res_Adult %>% select(-delabel) %>% arrange(desc(input_list_raw_log10p)) %>% head(5) %>%
    mutate(delabel = Tissue_cell_type_name)%>% select(Tissue_cell_type_name, delabel)
res_Adult = res_Adult  %>% select(-delabel)%>% left_join(label8)

label9 = res_Fetal %>% select(-delabel) %>% arrange(desc(input_list_raw_log10p)) %>% head(5) %>%
    mutate(delabel = Tissue_cell_type_name)%>% select(Tissue_cell_type_name, delabel)
res_Fetal = res_Fetal  %>% select(-delabel)%>% left_join(label9)


p8_clean <- ggplot(res_Adult, aes(x=Organ_system, y=input_list_raw_log10p,col=Organ_system)) + ggtitle("Adult Organ system") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res_Adult$input_list_raw_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_raw_pvalue_clean_Adult_Organ_system.pdf"),12,8)
print(p8)
dev.off()

p8 <- ggplot(res_Adult, aes(x=Organ_system, y=input_list_raw_log10p,col=Organ_system, label = delabel)) + geom_text_repel(size = 3)+ ggtitle("Adult Organ system") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res_Adult$input_list_raw_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_raw_pvalue_Adult_Organ_system.pdf"),12,8)
print(p8)
dev.off()

p9_clean <- ggplot(res_Fetal, aes(x=Organ_system, y=input_list_raw_log10p,col=Organ_system)) + ggtitle("Fetal Organ system") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res_Fetal$input_list_raw_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_raw_pvalue_clean_Fetal_Organ_system.pdf"),12,8)
print(p9)
dev.off()

p9 <- ggplot(res_Fetal, aes(x=Organ_system, y=input_list_raw_log10p,col=Organ_system, label = delabel))+ geom_text_repel(size = 3) + ggtitle("Fetal Organ system") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res_Fetal$input_list_raw_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_raw_pvalue_Fetal_Organ_system.pdf"),12,8)
print(p9)
dev.off()

##top 20 general cell types##

res_ordered <- res[order(res$input_list_raw_log10p,decreasing =T),]
top_20_cell_type <- unique(res_ordered$General_cell_type)[1:20]
cell_type_idx <- {}
for (i in 1:length(top_20_cell_type)){
	cell_type_idx <- c(cell_type_idx, which(res_ordered$General_cell_type==top_20_cell_type[i]))
}
res_ordered_top20_CT <- res_ordered[cell_type_idx,]
res_ordered_top20_CT$General_cell_type <- ordered(res_ordered_top20_CT$General_cell_type,levels = unique(res_ordered_top20_CT[order(res_ordered_top20_CT$input_list_raw_log10p,decreasing=T),6]))
label10 = res_ordered_top20_CT  %>% select(-delabel)%>% arrange(desc(input_list_raw_log10p)) %>% head(5) %>%
    mutate(delabel = Specific_cell_type)%>% select(Tissue_cell_type_name, delabel)
res_ordered_top20_CT = res_ordered_top20_CT  %>% select(-delabel)%>% left_join(label10)

p10_clean <- ggplot(res_ordered_top20_CT, aes(x=General_cell_type, y=input_list_raw_log10p,col=General_cell_type)) + ggtitle("Top20 general cell types") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res_ordered_top20_CT$input_list_raw_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_raw_pvalue_clean_Top20_general_cell_type.pdf"),12,8)
print(p10_clean)
dev.off()

p10 <- ggplot(res_ordered_top20_CT, aes(x=General_cell_type, y=input_list_raw_log10p,col=General_cell_type, label = delabel)) +geom_text_repel(size = 3) + ggtitle("Top20 general cell types") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(res_ordered_top20_CT$input_list_raw_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_raw_pvalue_Top20_general_cell_type.pdf"),12,8)
print(p10)
dev.off()

###########
##Finish!##
###########
