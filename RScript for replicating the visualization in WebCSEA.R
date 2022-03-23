##Rscript for replicating the visualization in WebCSEA (https://bioinfo.uth.edu/webcsea/).
##Here we use the Example gene list in the WebCSEA from Hepatocyte signature genes
library(ggplot2)
input_path <- "/PATH/TO/WebCSEA/DOWNLOAD/RESULTS/"
output_path <- "/PATH/TO/WORK/DIRECTORY/"
setwd(output_path)
##please find the WebCSEA.all_pvalue.txt file from the Download results zip folder.
heatmap_combine <- read.delim(paste0(output_path,"WebCSEA.all_pvalue.txt"),as.is=T,check.names=F)
##please find the curated organ,tissue,cell type information from the Download results zip folder.
organ_ref <- read.delim("organ_system_tissue_cell_type.txt",as.is=T,check.names=F)
##match the order of reference TCs and CSEA results, since some cells without CSEA results have been filtered out.
idx <- match(colnames(heatmap_combine),organ_ref$Tissue_cell_type_name)
if ( sum(is.na(idx))>=1){
idx <- idx[!is.na(idx)]
}else{
}

organ_ref_with_p <- data.frame(organ_ref[idx,],t(heatmap_combine),check.names=F)
organ_ref_with_p_logp <- data.frame(organ_ref_with_p,-log10( unlist(heatmap_combine[4,])),check.names=F)
colnames(organ_ref_with_p_logp)[ncol(organ_ref_with_p_logp)] <- "input_list_combined_log10p"
colnames(organ_ref_with_p_logp)[4] <- "Tissue"
p1 <- ggplot(organ_ref_with_p_logp, aes(x=Organ_system, y=input_list_combined_log10p,col=Organ_system)) + ggtitle("Organ system") + 
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(organ_ref_with_p_logp$input_list_combined_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"Organ_system.pdf"),12,8)
print(p1)
dev.off()

##top20 tissue cell_types result##
organ_ref_with_p_logp_ordered <- organ_ref_with_p_logp[order(organ_ref_with_p_logp$input_list_combined_log10p,decreasing =T),]
top_20_tissue <- unique(organ_ref_with_p_logp_ordered$Tissue)[1:20]
tissue_idx <- {}
for (i in 1:length(top_20_tissue)){
	tissue_idx <- c(tissue_idx, which(organ_ref_with_p_logp_ordered$Tissue==top_20_tissue[i]))
}
organ_ref_with_p_logp_top20 <- organ_ref_with_p_logp_ordered[tissue_idx,]
organ_ref_with_p_logp_top20$Tissue <- ordered(organ_ref_with_p_logp_top20$Tissue,levels = unique(organ_ref_with_p_logp_top20[order(organ_ref_with_p_logp_top20$input_list_combined_log10p,decreasing=T),4]))
p2 <- ggplot(organ_ref_with_p_logp_top20, aes(x=Tissue, y=input_list_combined_log10p,col=Tissue)) + ggtitle("Top20 tissues") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(organ_ref_with_p_logp_top20$input_list_combined_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"Top20_tissues.pdf"),12,8)
print(p2)
dev.off()
##Adult and Fetal organ_ref_with_p_logp
organ_ref_with_p_logp_Adult <- organ_ref_with_p_logp[organ_ref_with_p_logp$Developmental_stage=="Adult",]
organ_ref_with_p_logp_Fetal <- organ_ref_with_p_logp[organ_ref_with_p_logp$Developmental_stage=="Fetal",]

p3 <- ggplot(organ_ref_with_p_logp_Adult, aes(x=Organ_system, y=input_list_combined_log10p,col=Organ_system)) + ggtitle("Adult Organ system") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(organ_ref_with_p_logp_Adult$input_list_combined_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"Adult_Organ_system.pdf"),12,8)
print(p3)
dev.off()

p4 <- ggplot(organ_ref_with_p_logp_Fetal, aes(x=Organ_system, y=input_list_combined_log10p,col=Organ_system)) + ggtitle("Fetal Organ system") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(organ_ref_with_p_logp_Fetal$input_list_combined_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"Fetal_Organ_system.pdf"),12,8)
print(p4)
dev.off()

##top 20 general cell types##
organ_ref_with_p_logp_ordered <- organ_ref_with_p_logp[order(organ_ref_with_p_logp$input_list_combined_log10p,decreasing =T),]
top_20_cell_type <- unique(organ_ref_with_p_logp_ordered$General_cell_type)[1:20]
cell_type_idx <- {}
for (i in 1:length(top_20_cell_type)){
	cell_type_idx <- c(cell_type_idx, which(organ_ref_with_p_logp_ordered$General_cell_type==top_20_cell_type[i]))
}
organ_ref_with_p_logp_top20_CT <- organ_ref_with_p_logp_ordered[cell_type_idx,]
organ_ref_with_p_logp_top20_CT$General_cell_type <- ordered(organ_ref_with_p_logp_top20_CT$General_cell_type,levels = unique(organ_ref_with_p_logp_top20_CT[order(organ_ref_with_p_logp_top20_CT$input_list_combined_log10p,decreasing=T),6]))
p5 <- ggplot(organ_ref_with_p_logp_top20_CT, aes(x=General_cell_type, y=input_list_combined_log10p,col=General_cell_type)) + ggtitle("Top20 general cell types") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(organ_ref_with_p_logp_top20_CT$input_list_combined_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"Top20_general_cell_type.pdf"),12,8)
print(p5)
dev.off()

##
##generate the raw pvalue plots
organ_ref <- read.delim("organ_system_tissue_cell_type.txt",as.is=T,check.names=F)
idx <- match(colnames(heatmap_combine),organ_ref$Tissue_cell_type_name)
if ( sum(is.na(idx))>=1){
idx <- idx[!is.na(idx)]
}else{
}

organ_ref_with_p_raw <- data.frame(organ_ref[idx,],t(heatmap_combine),check.names=F)
organ_ref_with_p_logp <- data.frame(organ_ref_with_p_raw,-log10( unlist(heatmap_combine[1,])),check.names=F)
colnames(organ_ref_with_p_logp)[ncol(organ_ref_with_p_logp)] <- "input_list_raw_log10p"
colnames(organ_ref_with_p_logp)[4] <- "Tissue"
p6 <- ggplot(organ_ref_with_p_logp, aes(x=Organ_system, y=input_list_raw_log10p,col=Organ_system)) + ggtitle("Organ system") + 
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(organ_ref_with_p_logp$input_list_raw_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_raw_pvalue_Organ_system.pdf"),12,8)
print(p6)
dev.off()

##top20 tissue cell_types result##
organ_ref_with_p_logp_ordered <- organ_ref_with_p_logp[order(organ_ref_with_p_logp$input_list_raw_log10p,decreasing =T),]
top_20_tissue <- unique(organ_ref_with_p_logp_ordered$Tissue)[1:20]
tissue_idx <- {}
for (i in 1:length(top_20_tissue)){
	tissue_idx <- c(tissue_idx, which(organ_ref_with_p_logp_ordered$Tissue==top_20_tissue[i]))
}
organ_ref_with_p_logp_top20 <- organ_ref_with_p_logp_ordered[tissue_idx,]
organ_ref_with_p_logp_top20$Tissue <- ordered(organ_ref_with_p_logp_top20$Tissue,levels = unique(organ_ref_with_p_logp_top20[order(organ_ref_with_p_logp_top20$input_list_raw_log10p,decreasing=T),4]))
p7 <- ggplot(organ_ref_with_p_logp_top20, aes(x=Tissue, y=input_list_raw_log10p,col=Tissue)) + ggtitle("Top20 tissues") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(organ_ref_with_p_logp_top20$input_list_raw_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_raw_pvalue_Top20_tissues.pdf"),12,8)
print(p7)
dev.off()
##Adult and Fetal organ_ref_with_p_logp
organ_ref_with_p_logp_Adult <- organ_ref_with_p_logp[organ_ref_with_p_logp$Developmental_stage=="Adult",]
organ_ref_with_p_logp_Fetal <- organ_ref_with_p_logp[organ_ref_with_p_logp$Developmental_stage=="Fetal",]

p8 <- ggplot(organ_ref_with_p_logp_Adult, aes(x=Organ_system, y=input_list_raw_log10p,col=Organ_system)) + ggtitle("Adult Organ system") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(organ_ref_with_p_logp_Adult$input_list_raw_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_raw_pvalue_Adult_Organ_system.pdf"),12,8)
print(p8)
dev.off()

p9 <- ggplot(organ_ref_with_p_logp_Fetal, aes(x=Organ_system, y=input_list_raw_log10p,col=Organ_system)) + ggtitle("Fetal Organ system") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(organ_ref_with_p_logp_Fetal$input_list_raw_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_raw_pvalue_Fetal_Organ_system.pdf"),12,8)
print(p9)
dev.off()

##top 20 general cell types##
organ_ref_with_p_logp_ordered <- organ_ref_with_p_logp[order(organ_ref_with_p_logp$input_list_raw_log10p,decreasing =T),]
top_20_cell_type <- unique(organ_ref_with_p_logp_ordered$General_cell_type)[1:20]
cell_type_idx <- {}
for (i in 1:length(top_20_cell_type)){
	cell_type_idx <- c(cell_type_idx, which(organ_ref_with_p_logp_ordered$General_cell_type==top_20_cell_type[i]))
}
organ_ref_with_p_logp_top20_CT <- organ_ref_with_p_logp_ordered[cell_type_idx,]
organ_ref_with_p_logp_top20_CT$General_cell_type <- ordered(organ_ref_with_p_logp_top20_CT$General_cell_type,levels = unique(organ_ref_with_p_logp_top20_CT[order(organ_ref_with_p_logp_top20_CT$input_list_raw_log10p,decreasing=T),6]))
p10 <- ggplot(organ_ref_with_p_logp_top20_CT, aes(x=General_cell_type, y=input_list_raw_log10p,col=General_cell_type)) + ggtitle("Top20 general cell types") +
    geom_jitter() + theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = margin(0.5,0.5,1,2, "cm")) + geom_hline(yintercept=-log10(0.05/1355), linetype="dashed", color = "red", size=1) + geom_hline(yintercept=-log10(0.001), color = "grey", size=1) + ylim(0, max(organ_ref_with_p_logp_top20_CT$input_list_raw_log10p) +1.2)
pdf(paste0(output_path,"/",run_name,"_raw_pvalue_Top20_general_cell_type.pdf"),12,8)
print(p10)
dev.off()


###########
##Finish!##
###########
