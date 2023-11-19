library(CellChat)
library(ggplot2)
library(ggalluvial)
library(Seurat)
library(svglite)
library(patchwork)
source("./cellchat/cell_function.R")
# 1. GSE149614 ---------------------------------------------------------------
## 1.1 load data ---------------------------------------------------------------
Tumor <- readRDS("./outdata/GSE149614/Tumor.rds")
Lymphocyte <- readRDS("./outdata/GSE149614/lymphocyte/Lymphocyte4.rds")
Myeloid <- readRDS("./outdata/GSE149614/myeloid/Myeloid4.rds")

Tumor$main <- "Cancer"
Tumor$fine <- "Cancer"
merge_data <- merge(Tumor, y = c(Lymphocyte,Myeloid),merge.data = TRUE)
merge_data<-NormalizeData(merge_data)
saveRDS(merge_data, file = "./cellchat/GSE149614/merge_data.rds")
## 1.2 cellchat_object ---------------------------------------------------------------
##提取高低风险组
hep_imm_low <- subset(merge_data, orig.ident %in% c("HCC01T","HCC02T", "HCC03T","HCC04T","HCC05T","HCC06T"))
hep_imm_high <- subset(merge_data, orig.ident %in% c("HCC07T","HCC08T", "HCC09T","HCC10T"))
save(hep_imm_low,file = "./cellchat/GSE149614/hep_imm_low.Rdata")
save(hep_imm_high,file = "./cellchat/GSE149614/hep_imm_high.Rdata")
##每个单独跑cellchat_object，构建cellchat_object对象
cellchat_l<-cellchat_object(hep_imm_low)
cellchat_h<-cellchat_object(hep_imm_high)
save(cellchat_l,file = "./cellchat/GSE149614/cellchat_l.Rdata")
save(cellchat_h,file = "./cellchat/GSE149614/cellchat_h.Rdata")
object.list <- list(C1 = cellchat_h,C2 = cellchat_l)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
save(cellchat,file = "./cellchat/GSE149614/cellchat.Rdata")
save(object.list,file = "./cellchat/GSE149614/object_list.Rdata")  

####
library(ggsci)
## 1.3 Source and target ---------------------------------------------------------------
###Source and target
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf("./plot/cellchat/incoming_outgoing_GSE149614.pdf",width=8,height=4)
patchwork::wrap_plots(plots = gg)+plot_layout(guides = "collect")
dev.off()

## 1.4 Up or down-regulated receptor ligand ---------------------------------------------------------------

##
##DE
source("./cellchat/cellchat_mod.R")
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "C1"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes_median(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0,
              thresh.p = 1)
net <- netMappingDEG_median(cellchat, features.name = features.name,thresh = 0.01)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- net[net$datasets=="C1" &  net$ligand.median>=0 & net$receptor.median>=0 & (net$ligand.pvalues<0.05 | net$receptor.pvalues<0.05),]
net.down <- net[net$datasets=="C2" &  net$ligand.median<=0 & net$receptor.median<=0 & (net$ligand.pvalues<0.05 | net$receptor.pvalues<0.05),]
net.up <- net[net$datasets=="C1" &  ((net$ligand.median>=0 & net$ligand.pvalues<0.05) | (net$receptor.median>=0 & net$receptor.pvalues<0.05)) & 
                !((net$ligand.median>=0 & net$ligand.pvalues<0.05) & (net$receptor.median<=0 & net$receptor.pvalues<0.05)) & 
                !((net$ligand.median<=0 & net$ligand.pvalues<0.05) & (net$receptor.median>=0 & net$receptor.pvalues<0.05)) ,]

net.down <-net[net$datasets=="C2" &  ((net$ligand.median<=0 & net$ligand.pvalues<0.05) | (net$receptor.median<=0 & net$receptor.pvalues<0.05)) & 
                 !((net$ligand.median>=0 & net$ligand.pvalues<0.05) & (net$receptor.median<=0 & net$receptor.pvalues<0.05)) & 
                 !((net$ligand.median<=0 & net$ligand.pvalues<0.05) & (net$receptor.median>=0 & net$receptor.pvalues<0.05)),]


gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat) 
pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

## 1.5 Cancer ------------------------------------------------------------------
library(ggh4x)
cell=2
#from cancer (C1)
a123<-plot_point(cellchat,pairLR.use=pairLR.use.up,sources.use=cell,
                 targets.use=c(1:10),comparison=c(1,2),name="C1_deg_up",
                 max.dataset=1,net.up=net.up)
cancer_C1 <- df_pre(a123$data)
p1 <- LR_plot(rbind(cancer_C1))

#from cancer (C2)
a123<-plot_point(cellchat,pairLR.use=pairLR.use.down,sources.use=cell,
                 targets.use=c(1:10),comparison=c(1,2),name="C2_deg_pro",
                 max.dataset=2,net.up=net.down)  
cancer_C2 <- df_pre(a123$data)
p2 <- LR_plot(rbind(cancer_C2))

#to cancer (C1)
a123<-plot_point(cellchat,pairLR.use=pairLR.use.up,sources.use=c(1:10),
                 targets.use=cell,comparison=c(1,2),name="C1_deg_pro",
                 max.dataset=1,net.up=net.up)
cancer_C1_to <- df_pre(a123$data,reverse = T)
p3 <- LR_plot(rbind(cancer_C1,cancer_C1_to))

#to cancer (C2)
a123<-plot_point(cellchat,pairLR.use=pairLR.use.down,sources.use=c(1:10),
                 targets.use=cell,comparison=c(1,2),name="C2_deg_pro",
                 max.dataset=2,net.up=net.down)
cancer_C2_to <- df_pre(a123$data,reverse = T)
p4 <- LR_plot(rbind(cancer_C2,cancer_C2_to))

pdf("./plot/cellchat/Cancer_C1_GSE149614.pdf",width=10,height=3.4)
p3
dev.off()

pdf("./plot/cellchat/Cancer_C2_GSE149614.pdf",width=7,height=3)
p4
dev.off()

## 1.6 Macrophage ------------------------------------------------------------------
cell=6
a123<-plot_point(cellchat,pairLR.use=pairLR.use.up,sources.use=cell,
                 targets.use=c(1:10),comparison=c(1,2),name="C1_deg_up",
                 max.dataset=1,net.up=net.up)
cancer_C1 <- df_pre(a123$data)
p1 <- LR_plot(rbind(cancer_C1))

a123<-plot_point(cellchat,pairLR.use=pairLR.use.down,sources.use=cell,
                 targets.use=c(1:10),comparison=c(1,2),name="C2_deg_pro",
                 max.dataset=2,net.up=net.down)  
cancer_C2 <- df_pre(a123$data)
p2 <- LR_plot(rbind(cancer_C1))

a123<-plot_point(cellchat,pairLR.use=pairLR.use.up,sources.use=c(1:10),
                 targets.use=cell,comparison=c(1,2),name="C1_deg_pro",
                 max.dataset=1,net.up=net.up)
cancer_C1_to <- df_pre(a123$data,reverse = T)
p3 <- LR_plot(rbind(cancer_C1,cancer_C1_to))

a123<-plot_point(cellchat,pairLR.use=pairLR.use.down,sources.use=c(1:10),
                 targets.use=cell,comparison=c(1,2),name="C2_deg_pro",
                 max.dataset=2,net.up=net.down)
cancer_C2_to <- df_pre(a123$data,reverse = T)
p4 <- LR_plot(rbind(cancer_C2,cancer_C2_to))

pdf("./plot/cellchat/Macrophage_C1_GSE149614.pdf",width=9,height=3.2)
p3
dev.off()

pdf("./plot/cellchat/Macrophage_C2_GSE149614.pdf",width=10,height=5)
p4
dev.off()

## 1.7 SPP1 circle plot ------------------------------------------------------------------
pathways.show <- c("SPP1") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

pdf("./plot/cellchat/circle_SPP1_GSE149614.pdf",width=8,height=4)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


## 1.8 SPP1 Mac---------------------------------------------------------------
Tumor <- readRDS("./outdata/GSE149614/Tumor.rds")
Lymphocyte <- readRDS("./outdata/GSE149614/lymphocyte/Lymphocyte4.rds")
Myeloid <- readRDS("./outdata/GSE149614/myeloid/Myeloid4.rds")

Tumor$main <- "Cancer"
Tumor$fine <- "Cancer"
merge_data <- merge(Tumor, y = c(Lymphocyte,Myeloid),merge.data = TRUE)
merge_data<-NormalizeData(merge_data)
merge_data$main[merge_data$fine=="M1_Mø_SPP1"] <- "M1_Mø_SPP1"
merge_data$main[merge_data$main=="Macrophage"] <- "Other_Mø"
saveRDS(merge_data, file = "./cellchat/GSE149614/SPP1/merge_data.rds")
### 1.8.1 cellchat_object---------------------------------------------------------------
##提取高低风险组
hep_imm_low <- subset(merge_data, orig.ident %in% c("HCC01T","HCC02T", "HCC03T","HCC04T","HCC05T","HCC06T"))
hep_imm_high <- subset(merge_data, orig.ident %in% c("HCC07T","HCC08T", "HCC09T","HCC10T"))
save(hep_imm_low,file = "./cellchat/GSE149614/SPP1/hep_imm_low.Rdata")
save(hep_imm_high,file = "./cellchat/GSE149614/SPP1/hep_imm_high.Rdata")
##每个单独跑cellchat_object，构建cellchat_object对象
cellchat_l<-cellchat_object(hep_imm_low)
cellchat_h<-cellchat_object(hep_imm_high)
save(cellchat_l,file = "./cellchat/GSE149614/SPP1/cellchat_l.Rdata")
save(cellchat_h,file = "./cellchat/GSE149614/SPP1/cellchat_h.Rdata")
object.list <- list(C1 = cellchat_h,C2 = cellchat_l)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
save(cellchat,file = "./cellchat/GSE149614/SPP1/cellchat.Rdata")
save(object.list,file = "./cellchat/GSE149614/SPP1/object_list.Rdata")  

##上下调受体配体
##差异受体配体 
##差异情况
source("./cellchat/cellchat_mod.R")
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "C1"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes_median(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0,
                                              thresh.p = 1)
net <- netMappingDEG_median(cellchat, features.name = features.name,thresh = 0.01)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- net[net$datasets=="C1" &  net$ligand.median>=0 & net$receptor.median>=0 & (net$ligand.pvalues<0.05 | net$receptor.pvalues<0.05),]
net.down <- net[net$datasets=="C2" &  net$ligand.median<=0 & net$receptor.median<=0 & (net$ligand.pvalues<0.05 | net$receptor.pvalues<0.05),]
net.up <- net[net$datasets=="C1" &  ((net$ligand.median>=0 & net$ligand.pvalues<0.05) | (net$receptor.median>=0 & net$receptor.pvalues<0.05)) & 
                !((net$ligand.median>=0 & net$ligand.pvalues<0.05) & (net$receptor.median<=0 & net$receptor.pvalues<0.05)) & 
                !((net$ligand.median<=0 & net$ligand.pvalues<0.05) & (net$receptor.median>=0 & net$receptor.pvalues<0.05)) ,]

net.down <-net[net$datasets=="C2" &  ((net$ligand.median<=0 & net$ligand.pvalues<0.05) | (net$receptor.median<=0 & net$receptor.pvalues<0.05)) & 
                 !((net$ligand.median>=0 & net$ligand.pvalues<0.05) & (net$receptor.median<=0 & net$receptor.pvalues<0.05)) & 
                 !((net$ligand.median<=0 & net$ligand.pvalues<0.05) & (net$receptor.median>=0 & net$receptor.pvalues<0.05)),]


gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat) 
pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

###
cell=6
a123<-plot_point(cellchat,pairLR.use=pairLR.use.up,sources.use=cell,
                 targets.use=c(1:10),comparison=c(1,2),name="C1_deg_up",
                 max.dataset=1,net.up=net.up)
cancer_C1 <- df_pre(a123$data)


a123<-plot_point(cellchat,pairLR.use=pairLR.use.down,sources.use=cell,
                 targets.use=c(1:10),comparison=c(1,2),name="C2_deg_pro",
                 max.dataset=2,net.up=net.down)  
cancer_C2 <- df_pre(a123$data)


a123<-plot_point(cellchat,pairLR.use=pairLR.use.up,sources.use=c(1:10),
                 targets.use=cell,comparison=c(1,2),name="C1_deg_pro",
                 max.dataset=1,net.up=net.up)
cancer_C1_to <- df_pre(a123$data,reverse = T)
p3 <- LR_plot(rbind(cancer_C1,cancer_C1_to))

a123<-plot_point(cellchat,pairLR.use=pairLR.use.down,sources.use=c(1:10),
                 targets.use=cell,comparison=c(1,2),name="C2_deg_pro",
                 max.dataset=2,net.up=net.down)
cancer_C2_to <- df_pre(a123$data,reverse = T)
p4 <- LR_plot(rbind(cancer_C2,cancer_C2_to))

pdf("./plot/cellchat/MaSPP1_C1_GSE149614.pdf",width=9,height=3)
p3
dev.off()

pdf("./plot/cellchat/MaSPP1_C2_GSE149614.pdf",width=10,height=5)
p4
dev.off()


# 2. GSE151530 ------------------------------------------------------------
## 2.1 load data & cellchat object -----------------------------------------------------------
source("./cellchat/cell_function.R")
Tumor <- readRDS("E:/project/HCC_immune_scRNA/outdata/GSE151530/Tumor.rds")
Myeloid2 <- readRDS("E:/project/HCC_immune_scRNA/outdata/GSE151530/Myeloid2.rds")
Ly3 <- readRDS("E:/project/HCC_immune_scRNA/outdata/GSE151530/Ly3.rds")
Tumor$main <- "Cancer"
Tumor$fine <- "Cancer"
merge_data <- merge(Tumor, y = c(Ly3,Myeloid2),merge.data = TRUE)
merge_data<-NormalizeData(merge_data)
hep_imm_high<- subset(merge_data, Sample %in% c("H30","H37","H49b","H65","H70","H72","H74"))
hep_imm_low <- subset(merge_data, Sample %in% setdiff(unique(merge_data$Sample), 
                                                      c("H30","H37","H49b","H65","H70","H72","H74")))
cellchat_l<-cellchat_object(hep_imm_low)
cellchat_h<-cellchat_object(hep_imm_high)
object.list <- list(C1 = cellchat_h,C2 = cellchat_l)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

saveRDS(merge_data, file = "./cellchat/GSE151530/merge_data.rds")
save(cellchat_l,file = "./cellchat/GSE151530/cellchat_l.Rdata")
save(cellchat_h,file = "./cellchat/GSE151530/cellchat_h.Rdata")
save(cellchat,file = "./cellchat/GSE151530/cellchat.Rdata")
save(object.list,file = "./cellchat/GSE151530/object_list.Rdata")  

## 2.2 Source and target ---------------------------------------------------------------
###Source and target
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf("./plot/cellchat/incoming_outgoing_GSE151530.pdf",width=8,height=4)
patchwork::wrap_plots(plots = gg)+plot_layout(guides = "collect")
dev.off()

## 2.3 Up or down-regulated receptor ligand ---------------------------------------------------------------

##DE
source("./cellchat/cellchat_mod.R")
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "C1"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes_median(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0,
                                              thresh.p = 1)
net <- netMappingDEG_median(cellchat, features.name = features.name,thresh = 0.01)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- net[net$datasets=="C1" &  net$ligand.median>=0 & net$receptor.median>=0 & (net$ligand.pvalues<0.05 | net$receptor.pvalues<0.05),]
net.down <- net[net$datasets=="C2" &  net$ligand.median<=0 & net$receptor.median<=0 & (net$ligand.pvalues<0.05 | net$receptor.pvalues<0.05),]
net.up <- net[net$datasets=="C1" &  ((net$ligand.median>=0 & net$ligand.pvalues<0.05) | (net$receptor.median>=0 & net$receptor.pvalues<0.05)) & 
                !((net$ligand.median>=0 & net$ligand.pvalues<0.05) & (net$receptor.median<=0 & net$receptor.pvalues<0.05)) & 
                !((net$ligand.median<=0 & net$ligand.pvalues<0.05) & (net$receptor.median>=0 & net$receptor.pvalues<0.05)) ,]

net.down <-net[net$datasets=="C2" &  ((net$ligand.median<=0 & net$ligand.pvalues<0.05) | (net$receptor.median<=0 & net$receptor.pvalues<0.05)) & 
                 !((net$ligand.median>=0 & net$ligand.pvalues<0.05) & (net$receptor.median<=0 & net$receptor.pvalues<0.05)) & 
                 !((net$ligand.median<=0 & net$ligand.pvalues<0.05) & (net$receptor.median>=0 & net$receptor.pvalues<0.05)),]


gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat) 
pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

## 2.4 Cancer ------------------------------------------------------------------
cell=2
#from cancer (C1)
a123<-plot_point(cellchat,pairLR.use=pairLR.use.up,sources.use=cell,
                 targets.use=c(1:10),comparison=c(1,2),name="C1_deg_up",
                 max.dataset=1,net.up=net.up)
cancer_C1 <- df_pre(a123$data)

#from cancer (C2)
a123<-plot_point(cellchat,pairLR.use=pairLR.use.down,sources.use=cell,
                 targets.use=c(1:10),comparison=c(1,2),name="C2_deg_pro",
                 max.dataset=2,net.up=net.down)  
cancer_C2 <- df_pre(a123$data)

#to cancer (C1)
a123<-plot_point(cellchat,pairLR.use=pairLR.use.up,sources.use=c(1:10),
                 targets.use=cell,comparison=c(1,2),name="C1_deg_pro",
                 max.dataset=1,net.up=net.up)
cancer_C1_to <- df_pre(a123$data,reverse = T)
p3 <- LR_plot(rbind(cancer_C1,cancer_C1_to))

#to cancer (C2)
a123<-plot_point(cellchat,pairLR.use=pairLR.use.down,sources.use=c(1:10),
                 targets.use=cell,comparison=c(1,2),name="C2_deg_pro",
                 max.dataset=2,net.up=net.down)
cancer_C2_to <- df_pre(a123$data,reverse = T)
p4 <- LR_plot(rbind(cancer_C2,cancer_C2_to))

pdf("./plot/cellchat/Cancer_C1_GSE151530.pdf",width=9,height=4.2)
p3
dev.off()

pdf("./plot/cellchat/Cancer_C2_GSE151530.pdf",width=6,height=2)
p4
dev.off()

## 2.5 Macrophage ------------------------------------------------------------------
cell=6
a123<-plot_point(cellchat,pairLR.use=pairLR.use.up,sources.use=cell,
                 targets.use=c(1:10),comparison=c(1,2),name="C1_deg_up",
                 max.dataset=1,net.up=net.up)
cancer_C1 <- df_pre(a123$data)

a123<-plot_point(cellchat,pairLR.use=pairLR.use.down,sources.use=cell,
                 targets.use=c(1:10),comparison=c(1,2),name="C2_deg_pro",
                 max.dataset=2,net.up=net.down)  
cancer_C2 <- df_pre(a123$data)

a123<-plot_point(cellchat,pairLR.use=pairLR.use.up,sources.use=c(1:10),
                 targets.use=cell,comparison=c(1,2),name="C1_deg_pro",
                 max.dataset=1,net.up=net.up)
cancer_C1_to <- df_pre(a123$data,reverse = T)
p3 <- LR_plot(rbind(cancer_C1,cancer_C1_to))

a123<-plot_point(cellchat,pairLR.use=pairLR.use.down,sources.use=c(1:10),
                 targets.use=cell,comparison=c(1,2),name="C2_deg_pro",
                 max.dataset=2,net.up=net.down)
cancer_C2_to <- df_pre(a123$data,reverse = T)
p4 <- LR_plot(rbind(cancer_C2,cancer_C2_to))

pdf("./plot/cellchat/Macrophage_C1_GSE151530.pdf",width=10,height=6)
p3
dev.off()

pdf("./plot/cellchat/Macrophage_C2_GSE151530.pdf",width=7,height=2.5)
p4
dev.off()

## 2.6 SPP1 circle plot ------------------------------------------------------------------
pathways.show <- c("SPP1") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

pdf("./plot/cellchat/circle_SPP1_GSE151530.pdf",width=8,height=4)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()



## 2.7 SPP1 Mac ---------------------------------------------------------------
source("./cellchat/cell_function.R")
Tumor <- readRDS("E:/project/HCC_immune_scRNA/outdata/GSE151530/Tumor.rds")
Myeloid2 <- readRDS("E:/project/HCC_immune_scRNA/outdata/GSE151530/Myeloid2.rds")
Ly3 <- readRDS("E:/project/HCC_immune_scRNA/outdata/GSE151530/Ly3.rds")
Tumor$main <- "Cancer"
Tumor$fine <- "Cancer"
merge_data <- merge(Tumor, y = c(Ly3,Myeloid2),merge.data = TRUE)
merge_data<-NormalizeData(merge_data)
merge_data$main[merge_data$seurat_clusters=="3"] <- "Mø_SPP1"
merge_data$main[merge_data$main=="Macrophage"] <- "Other_Mø"

hep_imm_high<- subset(merge_data, Sample %in% c("H30","H37","H49b","H65","H70","H72","H74"))
hep_imm_low <- subset(merge_data, Sample %in% setdiff(unique(merge_data$Sample), 
                                                      c("H30","H37","H49b","H65","H70","H72","H74")))
cellchat_l<-cellchat_object(hep_imm_low)
cellchat_h<-cellchat_object(hep_imm_high)
object.list <- list(C1 = cellchat_h,C2 = cellchat_l)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

saveRDS(merge_data, file = "./cellchat/GSE151530/SPP1/merge_data.rds")
save(cellchat_l,file = "./cellchat/GSE151530/SPP1/cellchat_l.Rdata")
save(cellchat_h,file = "./cellchat/GSE151530/SPP1/cellchat_h.Rdata")
save(cellchat,file = "./cellchat/GSE151530/SPP1/cellchat.Rdata")
save(object.list,file = "./cellchat/GSE151530/SPP1/object_list.Rdata")  

##上下调受体配体
##差异受体配体 
##差异情况
source("./cellchat/cellchat_mod.R")
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "C1"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes_median(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0,
                                              thresh.p = 1)
net <- netMappingDEG_median(cellchat, features.name = features.name,thresh = 0.01)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- net[net$datasets=="C1" &  net$ligand.median>=0 & net$receptor.median>=0 & (net$ligand.pvalues<0.05 | net$receptor.pvalues<0.05),]
net.down <- net[net$datasets=="C2" &  net$ligand.median<=0 & net$receptor.median<=0 & (net$ligand.pvalues<0.05 | net$receptor.pvalues<0.05),]
net.up <- net[net$datasets=="C1" &  ((net$ligand.median>=0 & net$ligand.pvalues<0.05) | (net$receptor.median>=0 & net$receptor.pvalues<0.05)) & 
                !((net$ligand.median>=0 & net$ligand.pvalues<0.05) & (net$receptor.median<=0 & net$receptor.pvalues<0.05)) & 
                !((net$ligand.median<=0 & net$ligand.pvalues<0.05) & (net$receptor.median>=0 & net$receptor.pvalues<0.05)) ,]

net.down <-net[net$datasets=="C2" &  ((net$ligand.median<=0 & net$ligand.pvalues<0.05) | (net$receptor.median<=0 & net$receptor.pvalues<0.05)) & 
                 !((net$ligand.median>=0 & net$ligand.pvalues<0.05) & (net$receptor.median<=0 & net$receptor.pvalues<0.05)) & 
                 !((net$ligand.median<=0 & net$ligand.pvalues<0.05) & (net$receptor.median>=0 & net$receptor.pvalues<0.05)),]


gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat) 
pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

###
cell=7
a123<-plot_point(cellchat,pairLR.use=pairLR.use.up,sources.use=cell,
                 targets.use=c(1:10),comparison=c(1,2),name="C1_deg_up",
                 max.dataset=1,net.up=net.up)
cancer_C1 <- df_pre(a123$data)
LR_plot(rbind(cancer_C1))

a123<-plot_point(cellchat,pairLR.use=pairLR.use.down,sources.use=cell,
                 targets.use=c(1:10),comparison=c(1,2),name="C2_deg_pro",
                 max.dataset=2,net.up=net.down)  
cancer_C2 <- df_pre(a123$data)

a123<-plot_point(cellchat,pairLR.use=pairLR.use.up,sources.use=c(1:10),
                 targets.use=cell,comparison=c(1,2),name="C1_deg_pro",
                 max.dataset=1,net.up=net.up)
cancer_C1_to <- df_pre(a123$data,reverse = T)

p3 <- LR_plot(rbind(cancer_C1,cancer_C1_to))

a123<-plot_point(cellchat,pairLR.use=pairLR.use.down,sources.use=c(1:10),
                 targets.use=cell,comparison=c(1,2),name="C2_deg_pro",
                 max.dataset=2,net.up=net.down)
cancer_C2_to <- df_pre(a123$data,reverse = T)
p4 <- LR_plot(rbind(cancer_C2,cancer_C2_to))

pdf("./plot/cellchat/MaSPP1_C1_GSE151530.pdf",width=12,height=5)
p3
dev.off()

pdf("./plot/cellchat/MaSPP1_C2_GSE151530.pdf",width=15,height=4)
p4
dev.off()

# 3. SPP1 single-cell -----------------------------------------------------------------

## 3.1 GSE149614 -----------------------------------------------------------
##load data
merge_data <- readRDS("./cellchat/GSE149614/merge_data.rds")
unique(merge_data$fine)
merge_data$subtype<- ifelse(merge_data$Sample %in% c("HCC01T","HCC02T","HCC03T",
                                             "HCC04T","HCC05T","HCC06T"),"C2","C1")

source("./code/Function.R")
library(dittoSeq)
library(ggpubr)
p1<-dittoPlot(merge_data, "SPP1",group.by = "subtype",split.by = "fine",
              plots = c("boxplot"),
              boxplot.width = 0.7)+
  facet_wrap2(vars(fine),nrow = 3)+
  scale_fill_manual(values = c("C1"="#CD1818","C2"="#0F2C67"))+
  stat_compare_means(aes(group=grouping),vjust=0.7,
                     label = "p.signif", label.x.npc="middle",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1,1),
                                      symbols = c("***", "**", "*", "+","ns")))

plotdata<-p1$data
p2 <- plot.violin(plotdata,x=fine,y=var.data,color = grouping,
            group=grouping,add = c("boxplot","shadow"),width.b = 0.5,palette.col = c("#CD1818","#0F2C67"),
            shadow.col = c("grey40","white"),title = "SPP1")+
  theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1))
pdf("./plot/SPP1/SPP1_celltype_GSE149614.pdf",height = 3.5,width =10 )
p2
dev.off()

##
library(ggh4x)
Tumor <- subset(merge_data, fine=="Cancer")
plotdata <- data.frame(subtype=Tumor$subtype, ITGA5=Tumor@assays$RNA@data["ITGA5",], ITGB1=Tumor@assays$RNA@data["ITGB1",])
plotdata <- reshape2::melt(plotdata)
plotdata$Type <- "Cancer (GSE149614)"
p5<-plot.violin(plotdata,x=variable,y=value,color = subtype,
            group=subtype,add = c("boxplot","shadow"),width.b = 0.5,palette.col = c("#CD1818","#0F2C67"),
            shadow.col = c("grey40","white"),
            theme = theme_bw(base_size = 16),
            x_label="")+
  facet_wrap2("Type",scales ="free_y")  +
  theme(plot.margin = margin(t = 1,unit = 'cm'))+
  coord_cartesian(clip = 'off')
  

## 3.1 GSE151530 -----------------------------------------------------------

merge_data <- readRDS("./cellchat/GSE151530/merge_data.rds")
merge_data$seurat_clusters <- paste0(merge_data$main,"-",merge_data$seurat_clusters)
unique(paste0(merge_data$main,"-",merge_data$seurat_clusters))
merge_data$subtype<- ifelse(merge_data$Sample %in% c("H30","H37","H49b","H65","H70","H72","H74"),"C1","C2")
merge_data$seurat_clusters[is.na(merge_data$seurat_clusters)] <- "Cancer"

source("./code/Function.R")
library(dittoSeq)
library(ggh4x)
p1<-dittoPlot(merge_data, "SPP1",group.by = "subtype",split.by = "seurat_clusters",
              plots = c("boxplot"),
              boxplot.width = 0.7)+
  facet_wrap2(vars(seurat_clusters),nrow = 3)+
  scale_fill_manual(values = c("C1"="#CD1818","C2"="#0F2C67"))+
  stat_compare_means(aes(group=grouping),vjust=0.7,
                     label = "p.signif", label.x.npc="middle",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1,1),
                                      symbols = c("***", "**", "*", "+","ns")))

plotdata<-p1$data
p2 <- plot.violin(plotdata,x=seurat_clusters,y=var.data,color = grouping,
                  group=grouping,add = c("boxplot","shadow"),width.b = 0.5,palette.col = c("#CD1818","#0F2C67"),
                  shadow.col = c("grey40","white"),title = "SPP1")+
  theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1))
pdf("./plot/SPP1/SPP1_celltype_GSE151530.pdf",height = 3,width =9 )
p2
dev.off()


##
library(ggh4x)
Tumor <- subset(merge_data, seurat_clusters=="Cancer")
plotdata <- data.frame(subtype=Tumor$subtype, ITGA5=Tumor@assays$RNA@data["ITGA5",], ITGB1=Tumor@assays$RNA@data["ITGB1",])
plotdata <- reshape2::melt(plotdata)
plotdata$Type <- "Cancer (GSE151530)"
p6<-plot.violin(plotdata,x=variable,y=value,color = subtype,
                group=subtype,add = c("boxplot","shadow"),width.b = 0.5,palette.col = c("#CD1818","#0F2C67"),
                shadow.col = c("grey40","white"),
                theme = theme_bw(base_size = 16),
                x_label="")+
  facet_wrap2("Type",scales ="free_y")  +
  theme(plot.margin = margin(t = 1,unit = 'cm'))+
  coord_cartesian(clip = 'off')

library(patchwork)
pdf("./plot/SPP1/ITGA.pdf", width = 6, height = 3)
  p5+p6+plot_layout(guides = "collect")
dev.off()

