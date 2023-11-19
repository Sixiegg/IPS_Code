library(Seurat)
library(Matrix)
library(tidyverse)
library(patchwork)

# 1.Initial processing HCC_GSE149614 Lu_2022-----------------------------------------------------------
load("./initial_data/GSE149614/HCC.Rdata")
hcc.data <- CreateSeuratObject(counts = HCC,min.cells = 1)
hcc.data$Sample <- as.character(hcc.data$orig.ident)
hcc.data$Dataset <- "Lu_2022"
##按样本拆分
GSE149614 <- SplitObject(hcc.data,split.by='Sample')
saveRDS(GSE149614, file = "./data/GSE149614.rds")

# 2.DoubletFinder (remove double cell) -----------------------------------------------------------
library(DoubletFinder)
GSE149614 <- readRDS("./data/GSE149614.rds")
Doublet <- function(test, dim.usage=20){
  seu <- NormalizeData(test) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:dim.usage) %>%
    FindNeighbors(dims = 1:dim.usage) %>%
    FindClusters(resolution =0.5) %>%
    identity()
  
  ## pK Identification (no ground-truth)
  sweep.res.list <- paramSweep_v3(seu , PCs = 1:dim.usage, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  ##提取最优PK
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  ## Homotypic Doublet Proportion Estimate 
  annotations <- seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  DoubletRate = ncol(test)*8*1e-6
  nExp_poi <- round(DoubletRate*nrow(seu@meta.data))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies 
  seu <- doubletFinder_v3(seu, PCs = 1:dim.usage, pN = 0.25, 
                          pK = mpK, nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
  colnames(seu@meta.data)[ncol(seu@meta.data)] = "doublet_info"
  seu1 <- subset(seu, doublet_info=="Singlet")
  test1 <- CreateSeuratObject(counts = seu1@assays$RNA@counts,
                              meta.data = seu1@meta.data[,1:5],
                              min.cells = 1)
  return(test1)
}
test1 <- list(a=GSE149614$HCC01T)
GSE149614_dou_fil <- lapply(test1, function(x){Doublet(x)})


GSE149614_dou_fil <- lapply(GSE149614, function(x){Doublet(x)})
saveRDS(GSE149614_dou_fil, file = "./data/GSE149614_dou_fil.rds")

# 3.QC ----------------------------------------------------------------------
GSE149614_dou_fil <- readRDS("./data/GSE149614_dou_fil.rds")
merge_test <- merge(GSE149614_dou_fil[[1]], 
                    y=GSE149614_dou_fil[2:length(GSE149614_dou_fil)])

merge_test[["percent.mt"]] <- PercentageFeatureSet(merge_test, pattern = "^MT-")
merge_test[["nCount_RNA"]] = colSums(x = merge_test, slot = "counts")  # nCount_RNA
merge_test[["nFeature_RNA"]] = colSums(x = GetAssayData(object = merge_test, slot = "counts") > 0)

VlnPlot(object = merge_test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hcc.data1 <- subset(merge_test,percent.mt<=10 & nFeature_RNA >= 600 & nCount_RNA<=50000 & nFeature_RNA<=6000)
VlnPlot(object = hcc.data1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(hcc.data1$Sample)

saveRDS(hcc.data1, file = "./outdata/GSE149614/hcc.data1.rds")

# 4.First clustering ---------------------------------------------------------------------
library(Seurat)
library(harmony)
#merge_test1<-hcc.data1
##读入数据
merge_test1 <- readRDS("./outdata/GSE149614/hcc.data1.rds")
merge_test1 <- NormalizeData(merge_test1) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
  ScaleData()
merge_test1 <- RunPCA(merge_test1,verbose = F, npcs = 40)
merge_test1 <- merge_test1 %>% RunHarmony("orig.ident", plot_convergence = F,lambda=1)

merge_test1 <- merge_test1 %>% 
  # RunTSNE(reduction = "harmony", dims = 1:20) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20)

merge_test1 <- merge_test1 %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20)

# merge_test1 <- merge_test1 %>%
#   FindClusters(resolution = seq(from=0.5,by=.1,length=7)) %>%
#   identity()

merge_test1 <- merge_test1 %>% 
  FindClusters(resolution = 0.6) %>%
  identity()
saveRDS(merge_test1, file = "./outdata/GSE149614/merge_test1.rds")


##可视化聚类图
library(dittoSeq)
dittoDimPlot(merge_test1, reduction.use = "umap", var = "RNA_snn_res.0.6",
             do.label=T,size =0.2,opacity =0.8,legend.size=8)

FeaturePlot(merge_test1,reduction = "umap",features = c("CD3D","CD79A","TPSAB1","CD68","MZB1"))
FeaturePlot(merge_test1,reduction = "umap",features = c("CD3D","CD68","NKG7","CD79A"))
FeaturePlot(merge_test1,reduction = "umap",features = c("ALB","FGG","APOC3","FABP1"))
FeaturePlot(merge_test1,reduction = "umap",features = c("ALB","FGG","EPCAM","MZB1"))
FeaturePlot(merge_test1,reduction = "umap",features = c("ALB","FGG","EPCAM","MZB1"))
FeaturePlot(merge_test1,reduction = "umap",features = c("CD3D"))

##气泡图
# mac_name = list(
#   EPI_or_HEP = c("ALB","FGG","APOC3","FABP1","EPCAM"), #,"KRT19" 
#   Fibroblast =  c("COL1A2","ACTA2","COL3A1","MYH11"),
#   Endothelial = c("ENG","CDH5","VWF"),
#   Lymphocyte =  c("CD3D","NKG7","MZB1","CD79A"),
#   Myeloid =  c("CD68","C1QA","TPSAB1","CPA3")
# )
mac_name = list(
  EPI_or_HEP = c("ALB","FGG","APOC3","FABP1","EPCAM"), #,"KRT19"
  Fibroblast =  c("COL1A2","ACTA2","COL3A1","MYH11"),
  Endothelial = c("ENG","CDH5","VWF"),
  Lymphocyte =  c("CD3D","NKG7","MZB1","CD79A"),
  Myeloid =  c("CD68","C1QA","CD14")
)


# 生成矩阵气泡图，展示各类型巨噬细胞子集的表达情况
DotPlot(merge_test1, features = mac_name,
        assay='RNA' ,group.by = 'RNA_snn_res.0.6') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


merge_test2<-merge_test1
## 调整气泡图
# 重新排列聚类标签(Myeloid$dot)，并生成一个矩阵气泡图(Dot Plot)，展示不同类型巨噬细胞的表达谱及其在各个聚类中的表达情况
merge_test2$dot<-merge_test2$RNA_snn_res.0.6
merge_test2$dot<-factor(merge_test2$dot,levels=as.character(c(1,2,6,9,10,14,8,7,
                                                              3,4,11,12,13,15,
                                                              0,5,16,17)))
p2 <- DotPlot(merge_test2, features = mac_name,
            assay='RNA' ,group.by = 'dot') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p2
pdf(file = "./plot/main cluster/dotplot.pdf", width = 8,height = 5.5)
p2
dev.off()

##注释
# Main <- as.character(merge_test2$RNA_snn_res.0.8)
# Main[Main %in% c("1","2","6","9","10","14","15")]="Tumor"
# Main[Main %in% c("7")]="Fibroblast"
# Main[Main %in% c("8","19","21")]="Endothelial"
# Main[Main %in% c("3","4","11","12","13","16")]="Lymphocyte"
# Main[Main %in% c("0","5","17","18","20","22")]="Myeloid"
# merge_test2$Main <- Main

Main <- as.character(merge_test2$RNA_snn_res.0.6)
Main[Main %in% c("1","2","6","9","10","14")]="Tumor"
Main[Main %in% c("8")]="Fibroblast"
Main[Main %in% c("7")]="Endothelial"
Main[Main %in% c("3","4","11","12","13","15")]="Lymphocyte"
Main[Main %in% c("0","5","16","17")]="Myeloid"
merge_test2$Main <- Main


p1 <- dittoDimPlot(merge_test2, reduction.use = "umap", var = "RNA_snn_res.0.6",
                   do.label=T,size =0.2,opacity =0.8,legend.size=8)
p1
pdf(file = "./plot/main cluster/seurat cluster.pdf", width = 7,height = 6)
p1
dev.off()

p2 <- dittoDimPlot(merge_test2, reduction.use = "umap", var = "Main",
             do.label=T,size =0.2,opacity =0.8,legend.size=8)
p2
pdf(file = "./plot/main cluster/main cluster.pdf", width = 7,height = 6)
p2
dev.off()

##提取相关细胞群

Tumor <- subset(merge_test2, Main=="Tumor")
Fibroblast <- subset(merge_test2, Main=="Fibroblast")
Endothelial <- subset(merge_test2, Main=="Endothelial")
Lymphocyte <- subset(merge_test2, Main=="Lymphocyte")
Myeloid <- subset(merge_test2, Main=="Myeloid")

saveRDS(merge_test2, file = "./outdata/GSE149614/merge_test2.rds")
saveRDS(Tumor, file = "./outdata/GSE149614/Tumor.rds")
saveRDS(Fibroblast, file = "./outdata/GSE149614/Fibroblast.rds")
saveRDS(Endothelial, file = "./outdata/GSE149614/Endothelial.rds")
saveRDS(Lymphocyte, file = "./outdata/GSE149614/Lymphocyte.rds")
saveRDS(Myeloid, file = "./outdata/GSE149614/Myeloid.rds")


# 5. Sample subtype ------------------------------------------------------------
library(Seurat)
library(data.table)
source("./code/function.R")
library(GSVA)
library(ggpubr)
library(patchwork)
source("./code/validation.R")
library(dittoSeq)
result <- readRDS("./pair/result.rds")
result3 <- readRDS("./pair/result3.rds")
##sample subtype
{
  Tumor <- readRDS("./outdata/GSE149614/Tumor.rds")
  data<-Tumor_pair(Tumor,result,result3)
  data<-na.omit(data)
  subtype <- Tumor_class(Tumor,data,group.by = "Sample")
  saveRDS(subtype, file = "./pair/GSE149614_subtype.rds")
}

# 6. Lymphocyte ------------------------------------------------------------
## 6.1 Lymphocyte QC ----------------------------------------------------------
Lymphocyte <- readRDS("./outdata/GSE149614/Lymphocyte.rds")
Lymphocyte1 <- NormalizeData(Lymphocyte) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  #ScaleData(vars.to.regress=c("percent.mt","percent.ribo"))
  ScaleData() %>% 
  RunPCA(verbose = F, npcs = 40)
  
Lymphocyte1 <- Lymphocyte1 %>% RunHarmony("orig.ident", plot_convergence = F,
                                          max.iter.harmony=10)
  
Lymphocyte1 <- Lymphocyte1 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.7) %>% 
  identity()

MySeuratWrappers::VlnPlot(Lymphocyte1, features =c("CD3D","ALB","FGG","APOC3","FABP1","IRF7","TPSAB1"),
                          stacked=T,pt.size=0,cols =dittoColors(),
                          group.by = "RNA_snn_res.0.7")+NoLegend()

DimPlot(Lymphocyte1,reduction ="umap",label = T)

##11 dou，14 pDC,13 Mast
contain <- setdiff(as.character(0:14), c("11", "13", "14"))

Lymphocyte2 <- subset(Lymphocyte1, RNA_snn_res.0.7 %in% contain)
pDC <- subset(Lymphocyte1, RNA_snn_res.0.7 %in% as.character(c(14)))
mast <- subset(Lymphocyte1, RNA_snn_res.0.7 %in% as.character(c(13)))
Lymphocyte3 <- CreateSeuratObject(Lymphocyte2@assays$RNA@counts, meta.data = Lymphocyte2@meta.data[,1:6])
saveRDS(Lymphocyte3, file = "./outdata/GSE149614/lymphocyte/Lymphocyte_filter.rds")
saveRDS(pDC, file = "./outdata/GSE149614/myeloid/pDC.rds")
saveRDS(mast, file = "./outdata/GSE149614/myeloid/mast.rds")

## 6.2 Lymphocyte clustering ----------------------------------------------------------
Lymphocyte <- readRDS("./outdata/GSE149614/lymphocyte/Lymphocyte_filter.rds")
Lymphocyte1 <- NormalizeData(Lymphocyte) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
  #ScaleData(vars.to.regress=c("percent.mt","percent.ribo"))
  ScaleData() %>% 
  RunPCA(verbose = F, npcs = 50)

Lymphocyte1 <- Lymphocyte1 %>% RunHarmony("orig.ident", plot_convergence = F,
                                          max.iter.harmony=10)

Lymphocyte2 <- Lymphocyte1 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DimPlot(Lymphocyte2,reduction ="umap",label = T)
FeaturePlot(Lymphocyte2,reduction = "umap",features = c("KLRF1","CD4", "MZB1","CD8A"))
FeaturePlot(Lymphocyte2,reduction = "umap",features = c("NKG7"))
FeaturePlot(Lymphocyte2,reduction = "umap",features = c("CD4","NKG7","CD3D"))
dittoDimPlot(Lymphocyte2, reduction.use = "umap", var = "RNA_snn_res.0.5",
             do.label=T,size =0.2,opacity =0.8,legend.size=8)

FeaturePlot(Lymphocyte2,reduction = "umap",features = c("MKI67"))

marker_genes=c("CD3D","CD3E","CD4","CD8A","CD8B",##T
               "GZMA","GZMB","GZMH","GZMK","PRF1","GNLY",##effctor t
               "FOXP3","CTLA4","LAYN","TIGIT",##treg
               "IL7R","CCR7","TCF7","LEF1","SELL","S100A4",##tn/tm
               "NKG7","KLRF1","FGFBP2","KLRC1",##NK
               "CD79A","JSRP1","JCHAIN","IGHG1","MZB1",##B
               "HAVCR2","LAG3","CD274","PDCD1",###exhaustion
               "TOP2A","CDK1","MKI67",##proliferation
               "IFNG","CD69","CXCR6"##other
)

# 生成平均热图
library(scRNAtoolVis)
 AverageHeatmap(
  object = Lymphocyte2,                          # 分析对象
  markerGene = marker_genes,                 # 标记基因
  group.by = "RNA_snn_res.0.5",                # 按照聚类结果进行分组
  #cluster.order = as.character(c(0,2,3,4,11,13, 1,5,8,7,9,6,10,12)), # 聚类顺序
  cluster_columns = T,                   # 列不进行聚类
  #column_split = c(rep("CD4+T", 6), rep("CD8+T", 3), rep("NK", 2), rep("B", 3)), # 列分组
  row_split = c(rep("a", 5), rep("b", 6), rep("c", 4), rep("d", 6), rep("e", 4), rep("f", 5),
                rep("g", 4), rep("h", 3), rep("i", 3)), # 行分组
  border = TRUE                              # 边框
)
saveRDS(Lymphocyte2, file = "./outdata/GSE149614/lymphocyte/Lymphocyte2.rds")
 
## 6.3 Lymphocyte mix clusters(4,6) ------------------------------------------------------------------
##4 cluster
cluster4 <- subset(Lymphocyte2, RNA_snn_res.0.5 %in% c("4"))
#cluster4 <- CreateSeuratObject(cluster4@assays$RNA@counts)
DimPlot(cluster4,reduction ="umap",label = T)
FeaturePlot(cluster4,reduction = "umap",features = c("CD4", "NKG7","MZB1"))

cluster4 <- NormalizeData(cluster4) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>% 
  #ScaleData(vars.to.regress=c("percent.mt","percent.ribo"))
  ScaleData() %>% 
  RunPCA(verbose = F, npcs = 40)

cluster4 <- cluster4 %>% RunHarmony("orig.ident", plot_convergence = F,
                                          max.iter.harmony=10)

cluster4 <- cluster4 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()
DimPlot(cluster4,reduction ="umap",label = T)
FeaturePlot(cluster4,reduction = "umap",features = c("MZB1", "NKG7","CD4"))
FeaturePlot(cluster4,reduction = "umap",features = c("CD8A","CD4"))
Main <- as.character(cluster4$RNA_snn_res.0.8)
Main[Main %in% c("0","2")]="CD8+T"
Main[Main %in% c("1","3")]="CD4+T"
Main[Main %in% c("4","5")]="B"
cluster4$Main <- Main
DimPlot(cluster4,reduction ="umap",label = T,group.by="Main")


saveRDS(cluster4, file = "./outdata/GSE149614/lymphocyte/cluster4_mix.rds")

##修改注释
ano <- as.character(cluster4$RNA_snn_res.0.8)
ano[ano %in% c("1", "3")]="A"
ano[ano %in% c("0", "2")]="B"
ano[ano %in% c("4", "5")]="C"

ano1 <- as.character(Lymphocyte2$RNA_snn_res.0.5)
ano1[ano1 %in% c("4")]<- 
  paste0(ano1[ano1 %in% c("4")],
         ano[match(rownames(Lymphocyte2@meta.data)[ano1 %in% c("4")], rownames(cluster4@meta.data))])
  
##6 cluster
cluster6 <- subset(Lymphocyte2, RNA_snn_res.0.5 %in% c("6"))
DimPlot(cluster6,reduction ="umap",label = T)
FeaturePlot(cluster6,reduction = "umap",features = c("CD4", "NKG7"))

cluster6 <- NormalizeData(cluster6) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  #ScaleData(vars.to.regress=c("percent.mt","percent.ribo"))
  ScaleData() %>% 
  RunPCA(verbose = F, npcs = 50)

cluster6 <- cluster6 %>% RunHarmony("orig.ident", plot_convergence = F,
                                    max.iter.harmony=10)

cluster6 <- cluster6 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

DimPlot(cluster6,reduction ="umap",label = T)
FeaturePlot(cluster6,reduction = "umap",features = c("CD8A", "NKG7","CD4"))
FeaturePlot(cluster6,reduction = "umap",features = c("CD8A","CD4"))
saveRDS(cluster6, file = "./outdata/GSE149614/lymphocyte/cluster6_mix.rds")

##修改注释
ano6 <- as.character(cluster6$RNA_snn_res.0.8)
ano6[ano6 %in% c("0")]="B"
ano6[ano6 %in% c("1", "2", "3")]="A"

#ano1 <- as.character(Lymphocyte2$RNA_snn_res.0.5)
ano1[ano1 %in% c("6")]<- 
  paste0(ano1[ano1 %in% c("6")],
         ano6[match(rownames(Lymphocyte2@meta.data)[ano1 %in% c("6")], rownames(cluster6@meta.data))])

ano_cluster<-ano1
ano_cluster <- factor(ano_cluster,levels = c("0","1","2","3","4A","4B","4C","5","6A","6B","7","8","9"))
Lymphocyte3 <- Lymphocyte2
Lymphocyte3$ano_cluster <- ano_cluster

DimPlot(Lymphocyte3,reduction ="umap",label = T)
dittoDimPlot(Lymphocyte3, reduction.use = "umap", var = "ano_cluster",
             do.label=T,size =0.2,opacity =0.8,legend.size=8)
saveRDS(Lymphocyte3, file = "./outdata/GSE149614/lymphocyte/Lymphocyte3.rds")

## 6.4 Lymphocyte plot ------------------------------------------------------------------
### 6.4.1 Lymphocyte marker plot ------------------------------------------------------------------
Lymphocyte3 <- readRDS("./outdata/GSE149614/lymphocyte/Lymphocyte3.rds")  
Idents(Lymphocyte3) <- Lymphocyte3$ano_cluster
#平均热图(论文出图版)
{
  cluster_type <- data.frame(cluster = c("L0","L1","L4A","L6A","L2","L4B","L6B","L8","L9","L3",
                               "L4C","L5","L7"), cell_type=c(rep("CD4+T",4), rep("CD8+T", 3), rep("NK",2),
                                                             rep("B", 4)))
  marker_genes=c("CD3D","CD3E","CD4","CD8A","CD8B",##T
                 "GZMA","GZMB","GZMH","GZMK","PRF1","GNLY",##effctor t
                 "FOXP3","CTLA4","LAYN","TNFRSF18","IL2RA",##treg
                 "IL7R","CCR7","TCF7","LEF1",##tn/tm
                 "NKG7","KLRF1","FGFBP2","KLRC1",##NK
                 "CD79A","JSRP1","JCHAIN","IGHG1","MZB1",##B
                 "HAVCR2","LAG3","CD274","PDCD1",###exhaustion
                 "TOP2A","CDK1","MKI67"##proliferation
  )
  markerdata <- data.frame(gene=marker_genes, genetype=c(rep("T",5),
                                                         rep("Effector T",6),
                                                         rep("Treg",5),
                                                         rep("Tn/Tm",4),
                                                         rep("NK",4),
                                                         rep("B",5),
                                                         rep("Exhaustion",4),
                                                         rep("Proliferation",3)))
  ##mean scale
  lymean <- AverageExpression(Lymphocyte3,features = marker_genes, group.by="ano_cluster")
  htdf <- t(scale(t(lymean$RNA), scale = T, center = T))
  colnames(htdf) <- paste0("L", colnames(htdf))
  htdf <- data.frame(gene=marker_genes,htdf)
  plotdata <- reshape2::melt(htdf)
  plotdata$genetype <- markerdata$genetype[match(plotdata$gene, markerdata$gene)]
  plotdata$celltype <- cluster_type$cell_type[match(plotdata$variable, cluster_type$cluster)]
  plotdata$celltype <- factor(plotdata$celltype, levels = c("CD4+T", "CD8+T", "NK", "B"))
  plotdata$gene <- factor(plotdata$gene, levels = marker_genes[36:1])
  plotdata$genetype <- factor(plotdata$genetype, levels = c("T", "Effector T", "Treg", "Tn/Tm",
                                                            "NK", "B", "Exhaustion", "Proliferation")[8:1])
  ##heatmap
  library(ggpubr)
  library(ggh4x)
  htCol = c("#0099CC", "white","#CC0033")
  col_fun <- colorRampPalette(htCol)(50)
  ly_heatmap <- ggplot(plotdata, aes(x = interaction(variable, celltype), y = interaction(gene, genetype), fill = value)) +
    geom_raster() +
    scale_fill_gradientn(colors = col_fun, guide = "colorbar", na.value = "white", limits = c(-2, 2), 
                         breaks = seq(-2, 2, 1), labels = seq(-2, 2, 1), 
                         oob = scales::squish, name = "value") +
    xlab(NULL) + ylab(NULL) +
    theme(panel.grid = element_blank(), panel.background = element_blank()) +
    geom_tile() +
    theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) +
    geom_vline(xintercept = c(4.5, 7.5, 9.5), size = 1.5, color = "white") +
    geom_hline(yintercept = c(36.5 - 5, 36.5 - 11, 36.5 - 16, 36.5 - 20, 36.5 - 24, 36.5 - 29, 36.5 - 33, 36.5 - 36),
               size = 1.5, color = "white") +
    scale_y_discrete(position = 'left') +
    scale_x_discrete(position = 'top') +
    guides(x = "axis_nested", y = "axis_nested") +
    font("xy", size = 11) +
    font("xy.text", size = 11, color = "black") +
    font("legend.text", size = 12) +
    font("legend.title", size = 12) +
    theme(axis.title.x = element_blank())
} 
  
### 6.4.2 Lymphocyte OR plot ------------------------------------------------------------------
##分类情况
  subtype <- readRDS("./pair/GSE149614_subtype.rds") 

##OR图
  ##计算OR值
{
  test.dist.table <- function(count.dist,min.rowSum=0)
  {
    count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
    sum.col <- colSums(count.dist)
    sum.row <- rowSums(count.dist)
    count.dist.tb <- as.data.frame(count.dist)
    setDT(count.dist.tb,keep.rownames=T)
    count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
    colnames(count.dist.melt.tb) <- c("rid","cid","count")
    count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
      this.row <- count.dist.melt.tb$rid[i]
      this.col <- count.dist.melt.tb$cid[i]
      this.c <- count.dist.melt.tb$count[i]
      other.col.c <- sum.col[this.col]-this.c
      this.m <- matrix(c(this.c,
                         sum.row[this.row]-this.c,
                         other.col.c,
                         sum(sum.col)-sum.row[this.row]-other.col.c),
                       ncol=2)
      res.test <- fisher.test(this.m)
      data.frame(rid=this.row,
                 cid=this.col,
                 p.value=res.test$p.value,
                 OR=res.test$estimate)
    }))
    count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                    by=c("rid","cid"))
    count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
    return(count.dist.melt.ext.tb)
  }
  
    metadata1 <- Lymphocyte3@meta.data
    temp <- metadata1[,c("orig.ident","ano_cluster"),drop=F]
    temp$ano_cluster <- paste0("L",temp$ano_cluster)
    temp$subtype <- subtype[match(temp$orig.ident, rownames(subtype)), 4]
  
    library(data.table)
    library(plyr)
    cell_type<-data.table(temp[,-1])
    a1<-cell_type$ano_cluster
    a2<-cell_type$subtype
    temp <- unclass(cell_type[,table(a1,a2)])
    temp <- temp[unique(a1),]
    
    ly<-test.dist.table(temp)
}
  
  ###OR plot
{
  ly1<-ly
  ly1$OR[ly1$OR>3]=3
  #ly1$cid<-factor(ly1$cid,levels = c("C2","C1"))
  ly1$dir <- ifelse(ly1$rid %in% ly1[ly1$cid=="C1",]$rid[ly1[ly1$cid=="C1",]$OR >1], 1, -1)
  ly1 <- data.frame(ly1,dirOR=ly1$OR *  ly1$dir)
  ly1 <- data.frame(ly1,plabel=ifelse(ly1$p.value<0.05, "P<0.05","P≥0.05"))
  ly1$plabel <- factor(ly1$plabel,levels = c("P≥0.05","P<0.05"))
  cluster_type <- data.frame(cluster = c("L0","L1","L4A","L6A","L2","L4B","L6B","L8","L9","L3",
                                         "L4C","L5","L7"), cell_type=c(rep("CD4+T",4), rep("CD8+T", 3), rep("NK",2),
                                                                       rep("B", 4)))
  ly1$celltype <- cluster_type$cell_type[match(ly1$rid, cluster_type$cluster)]
  ly1$celltype <- factor(ly1$celltype, levels = c("CD4+T", "CD8+T", "NK", "B"))
  
  ly_dot <- ggplot(ly1,aes(x = interaction(rid,celltype),y = cid,color = dirOR,size=plabel)) +
    geom_point(size=6) +
    #scale_size_continuous(range=c(6,1))+
    #scale_size(range = c(1, 6))+
    cowplot::theme_cowplot() +
    theme(axis.line = element_blank()) +
    #scale_colour_viridis_c( option = "D",alpha =1,begin=1,end = 0)+
    scale_colour_stepsn(
      #n.breaks=40,
      #guide="legend",
      breaks = c(-3,-2,-1,0,1,2,3),
      #breaks = c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75),
      #labels=c(-3,-2,-1,0,1,2,3),
      na.value = 'white',
      colors =c(
        colorRampPalette(colors = c("#0F2C67","#BAC2D3"))(3),
        colorRampPalette(colors = c("#F0BCBC","#CD1818"))(3)
      ),name="OR"
    )+
    guides(x="axis_nested")+
    font("xy.text", size = 11,color = "black")+
    theme(axis.ticks = element_blank(),axis.title.x=element_blank(),
          axis.text.x=element_blank(),axis.title.y=element_blank())
  
  
  library(aplot)
  pdf("./plot/lymphocyte/GSE149614/marker_heat_dot.pdf", height = 8,width = 6.3)
  ly_heatmap %>% insert_bottom(ly_dot,height = 0.08)
  dev.off()
}

## 6.5 marker gene Lymphocyte ano ------------------------------------------------------------------
marker<-FindAllMarkers(Lymphocyte3,only.pos = T,logfc.threshold = 0.5,min.pct = 0.2)
markers1 <- marker[!grepl("^RP[ SL ]", marker$gene, ignore.case = F),]
markers1 <- markers1[!grepl("^MT-", markers1$gene, ignore.case = F),]
saveRDS(markers1, file = "./outdata/GSE149614/lymphocyte/markers1.rds")

library(dplyr)
library(dittoSeq)
fine1 <- as.character(Lymphocyte3$ano_cluster)
fine<-recode(fine1,
             "0" = "L0_CD4_IL7R",
             "1" = "L1_Treg_IL2RA",
             "2" = "L2_CD8_GZMK",
             "3" = "L3_B_IGHG1",
             "4A" = "L4A_Treg_TOP2A",
             "4B" = "L4B_CD8_TOP2A",
             "4C" = "L4C_B_MYBL2",
             "5" = "L5_B_MZB1",
             "6A" = "L6A_Treg_TNFRSF18",
             "6B" = "L6B_CD8_LAG3",
             "7" = "L7_B_CD79A",
             "8" = "L8_NK_KLRC1",
             "9" = "L9_NK_FGFBP2"
)
main<-recode(fine1,
             "0" = "CD4",
             "1" = "Treg",
             "2" = "CD8",
             "3" = "B",
             "4A" = "Treg",
             "4B" = "CD8",
             "4C" = "B",
             "5" = "B",
             "6A" = "Treg",
             "6B" = "CD8",
             "7" = "B",
             "8" = "NK",
             "9" = "NK"
)
Lymphocyte4 <- Lymphocyte3
Lymphocyte4$fine<-fine
Lymphocyte4$main<-main
dittoDimPlot(Lymphocyte4, reduction.use = "umap", var = "main",
             do.label=F,size =0.2,opacity =0.8,legend.size=8)
dittoDimPlot(Lymphocyte4, reduction.use = "umap", var = "fine",
             do.label=F,size =0.2,opacity =0.8,legend.size=8)
saveRDS(Lymphocyte4, file = "./outdata/GSE149614/lymphocyte/Lymphocyte4.rds")

pdf("./plot/lymphocyte/GSE149614/ly_ano_plot.pdf", height = 4.5,width = 6)
dittoDimPlot(Lymphocyte4, reduction.use = "umap", var = "fine",
             do.label=F,size =0.2,opacity =0.8,legend.size=8)
dev.off()


# 7. Myeloid clustering ------------------------------------------------------------
## 7.1 merge data  ---------------------------------------------------------------------
library(Seurat)
library(Matrix)
library(tidyverse)
library(patchwork)
library(harmony)
Myeloid <- readRDS("./outdata/GSE149614/Myeloid.rds")
Myeloid <- CreateSeuratObject(Myeloid@assays$RNA@counts, meta.data = Myeloid@meta.data[,1:6])
pDC <- readRDS("./outdata/GSE149614/myeloid/pDC.rds")
pDC <- CreateSeuratObject(pDC@assays$RNA@counts, meta.data = pDC@meta.data[,1:6])
mast <- readRDS("./outdata/GSE149614/myeloid/mast.rds")
mast <- CreateSeuratObject(mast@assays$RNA@counts, meta.data = mast@meta.data[,1:6])

#
Myeloid1 <- merge(Myeloid, y = c(pDC, mast),merge.data = TRUE)
saveRDS(Myeloid1, file = "./outdata/GSE149614/myeloid/Myeloid1.rds")


## 7.2 clustering --------------------------------------------------------------------
Myeloid2 <- NormalizeData(Myeloid1) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  #ScaleData(vars.to.regress=c("percent.mt","percent.ribo"))
  ScaleData() %>% 
  RunPCA(verbose = F, npcs = 50)

Myeloid2 <- Myeloid2 %>% RunHarmony("orig.ident", plot_convergence = F,
                                    max.iter.harmony=10)

Myeloid3 <- Myeloid2 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:15) %>% 
  FindClusters(resolution = 0.7) %>% 
  identity()
DimPlot(Myeloid3,reduction ="umap",label = T)

saveRDS(Myeloid3, file = "./outdata/GSE149614/myeloid/Myeloid3.rds")

DimPlot(Myeloid3,reduction ="umap",label = T)
FeaturePlot(Myeloid3,reduction = "umap",features = c("CLEC10A","C1QA"))
FeaturePlot(Myeloid3,reduction = "umap",features = c("SPP1","C1QA"))
marker<-FindAllMarkers(Myeloid3,only.pos = T,min.pct = 0.1)
markers1 <- marker[!grepl("^RP[ SL ]", marker$gene, ignore.case = F),]
markers1 <- markers1[!grepl("^MT-", markers1$gene, ignore.case = F),]
#marker2 <- markers1[markers1$cluster=="1", ]

VlnPlot(object = Myeloid3, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


mac_name = list(
  Mac=c("LYZ","CD68","CD14","CD163","C1QA","C1QB","CXCL9"),
  Mono=c("VCAN","FCN1","S100A8","S100A9","FCGR3A"),
  Neutrophils =  c("FCGR3B","CXCR2","SLC25A37","G0S2","CXCR1","ADGRG3","PROK2","STEAP4","CMTM2" ),
  pDC = c("GZMB","CLEC4C","IRF7"),
  DC1 = c("CLEC9A","IRF8","IDO1","HLA-DPB1","HLA-DQB1"),
  DC2=c( "CD1C","CD1E","CD1A","CLEC10A"),
  Mast=c("TPSAB1","CPA3"),
  KC=c("CD5L","VCAM1","MARCO","TIMD4"),
  Other=c("SPP1","MMP9","CD3D","MKI67","TOP2A","PMEL","MLANA","MRC1")
)

mac_name = list(
  Mac=c("LYZ","CD68","CD14","CD163","C1QA","C1QB","CXCL9")
)
# 生成矩阵气泡图，展示各类型巨噬细胞子集的表达情况
DotPlot(Myeloid3, features = c("LYZ","CD68","CD14","CD163","C1QA","C1QB","CXCL9"),
        assay='RNA' ,group.by = 'seurat_clusters') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


## 7.3 marker plot --------------------------------------------------------------------
### 7.3.1 heatmap --------------------------------------------------------------------
##
Myeloid3 <- readRDS("./outdata/GSE149614/myeloid/Myeloid3.rds")
marker_genes=c("CD68","CD14","CD163","FCGR3A",##Macrophage
               "SPP1","MMP7","MMP9","ITGAM","VEGFA",##Tam
               "FCN1","LYZ","VCAN", 'S100A8', 'S100A9',##Monocyte
               'CXCL9','CXCL10',"CXCL11",##M1 like
               'MRC1',"TGFB1",'ARG1',"GPNMB",##M2 like
               "GZMB","IRF7","CLEC9A","IRF8","IDO1","CD1C","CD1E","CLEC10A",##DC
               "HLA-DRA","HLA-A","C1QA","C1QC",##APC&&CC
               "TPSAB1","CPA3",#mast
               "TOP2A","CDK1","MKI67",##Proliferation
               "APOC1","APOE","CCL5"##other
)

#平均热图(论文出图版)
{
  cluster_type <- data.frame(cluster =  paste0("M", as.character(c(0,2,4,8,1,6,3,5,7,9,11,10))), cell_type=c(rep("Mø", 6), rep("Mon", 2), rep("DC", 3), rep("Mast", 1)))
  marker_genes=c("CD68","CD14","CD163","FCGR3A",##Macrophage
                 "FCN1","LYZ","VCAN", 'S100A8', 'S100A9',##Monocyte
                 "GZMB","IRF7","CLEC9A","IRF8","IDO1","CD1C","CD1E","CLEC10A",##DC
                 "TPSAB1","CPA3",#mast
                 "SPP1","MMP7","MMP9","MMP12","SPARC",##EMR
                 "HLA-DRA","HLA-A","C1QA","C1QC",##APC&&CC
                 "TOP2A","CDK1","MKI67"##Proliferation
  )
  markerdata <- data.frame(gene=marker_genes, genetype=c(rep("Mø",4),
                                                         rep("Mon",5),
                                                         rep("DC",8),
                                                         rep("Mast",2),
                                                         rep("EMR",5),
                                                         rep("APC&&CC",4),
                                                         rep("Proliferation",3)))
  ##mean scale
  mymean <- AverageExpression(Myeloid3,features = marker_genes, group.by="seurat_clusters")
  htdf <- t(scale(t(mymean$RNA), scale = T, center = T))
  colnames(htdf) <- paste0("M", colnames(htdf))
  htdf <- data.frame(gene=marker_genes,htdf)
  plotdata <- reshape2::melt(htdf)
  plotdata$genetype <- markerdata$genetype[match(plotdata$gene, markerdata$gene)]
  plotdata$celltype <- cluster_type$cell_type[match(plotdata$variable, cluster_type$cluster)]
  plotdata$celltype <- factor(plotdata$celltype, levels = c("Mø", "Mon", "DC", "Mast"))
  plotdata$gene <- factor(plotdata$gene, levels = marker_genes[31:1])
  plotdata$genetype <- factor(plotdata$genetype, levels = c("Mø", "Mon", "DC", "Mast",
                                                            "EMR", "APC&&CC", "Proliferation")[7:1])
  plotdata$variable <- factor(plotdata$variable, levels = paste0("M", as.character(c(0,2,4,8,1,6,3,5,7,9,11,10))))
  
  ##heatmap
  library(ggpubr)
  library(ggh4x)
  htCol = c("#0099CC", "white","#CC0033")
  col_fun <- colorRampPalette(htCol)(50)
  my_heatmap <- 
    ggplot(plotdata, aes(x = interaction(variable, celltype), y = interaction(gene, genetype), fill = value)) +
    geom_raster() +
    scale_fill_gradientn(colors = col_fun, guide = "colorbar", na.value = "white", limits = c(-2, 2), 
                         breaks = seq(-2, 2, 1), labels = seq(-2, 2, 1), 
                         oob = scales::squish, name = "value") +
    xlab(NULL) + ylab(NULL) +
    theme(panel.grid = element_blank(), panel.background = element_blank()) +
    geom_tile() +
    theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) +
    geom_vline(xintercept = c(6.5, 8.5, 11.5,12.5), size = 1.5, color = "white") +
    geom_hline(yintercept = c(31.5 - 4, 31.5 - 9, 31.5 - 17, 31.5 - 19, 31.5 - 24, 31.5 - 28),
               size = 1.5, color = "white") +
    scale_y_discrete(position = 'left') +
    scale_x_discrete(position = 'top') +
    guides(x = "axis_nested", y = "axis_nested") +
    font("xy", size = 11) +
    font("xy.text", size = 11, color = "black") +
    font("legend.text", size = 12) +
    font("legend.title", size = 12) +
    theme(axis.title.x = element_blank())
} 
### 7.3.2 OR plot ------------------------------------------------------------------
##分类情况
subtype <- readRDS("./pair/GSE149614_subtype.rds") 

##OR图
##计算OR值
{
  test.dist.table <- function(count.dist,min.rowSum=0)
  {
    count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
    sum.col <- colSums(count.dist)
    sum.row <- rowSums(count.dist)
    count.dist.tb <- as.data.frame(count.dist)
    setDT(count.dist.tb,keep.rownames=T)
    count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
    colnames(count.dist.melt.tb) <- c("rid","cid","count")
    count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
      this.row <- count.dist.melt.tb$rid[i]
      this.col <- count.dist.melt.tb$cid[i]
      this.c <- count.dist.melt.tb$count[i]
      other.col.c <- sum.col[this.col]-this.c
      this.m <- matrix(c(this.c,
                         sum.row[this.row]-this.c,
                         other.col.c,
                         sum(sum.col)-sum.row[this.row]-other.col.c),
                       ncol=2)
      res.test <- fisher.test(this.m)
      data.frame(rid=this.row,
                 cid=this.col,
                 p.value=res.test$p.value,
                 OR=res.test$estimate)
    }))
    count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                    by=c("rid","cid"))
    count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
    return(count.dist.melt.ext.tb)
  }
  
  metadata1 <- Myeloid3@meta.data
  temp <- metadata1[,c("orig.ident","seurat_clusters"),drop=F]
  temp$seurat_clusters <- paste0("M",temp$seurat_clusters)
  temp$subtype <- subtype[match(temp$orig.ident, rownames(subtype)), 4]
  
  library(data.table)
  library(plyr)
  cell_type<-data.table(temp[,-1])
  a1<-cell_type$seurat_clusters
  a2<-cell_type$subtype
  temp <- unclass(cell_type[,table(a1,a2)])
  temp <- temp[unique(a1),]
  
  my<-test.dist.table(temp)
}

###OR plot
{
  my1<-my
  my1$OR[my1$OR>3]=3
  #my1$cid<-factor(my1$cid,levels = c("C2","C1"))
  my1$dir <- ifelse(my1$rid %in% my1[my1$cid=="C1",]$rid[my1[my1$cid=="C1",]$OR >1], 1, -1)
  my1 <- data.frame(my1,dirOR=my1$OR *  my1$dir)
  my1 <- data.frame(my1,plabel=ifelse(my1$p.value<0.05, "P<0.05","P≥0.05"))
  my1$plabel <- factor(my1$plabel,levels = c("P≥0.05","P<0.05"))
  cluster_type <- data.frame(cluster = c("M0","M2","M4","M8","M1","M6","M3","M5","M7","M9",
                                         "M11","M10"), cell_type=c(rep("Mø", 6), rep("Mon", 2), rep("DC", 3), rep("Mast", 1)))
  my1$celltype <- cluster_type$cell_type[match(my1$rid, cluster_type$cluster)]
  my1$celltype <- factor(my1$celltype, levels = c("Mø", "Mon", "DC", "Mast"))
  
  my_dot <- 
    ggplot(my1,aes(x = interaction(rid,celltype),y = cid,color = dirOR,size=plabel)) +
    geom_point() +
    scale_size_discrete(range=c(3,6))+
    #scale_size(range = c(1, 6))+
    cowplot::theme_cowplot() +
    theme(axis.line = element_blank()) +
    #scale_colour_viridis_c( option = "D",alpha =1,begin=1,end = 0)+
    scale_colour_stepsn(
      #n.breaks=40,
      #guide="legend",
      breaks = c(-3,-2,-1,0,1,2,3),
      #breaks = c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75),
      #labels=c(-3,-2,-1,0,1,2,3),
      na.value = 'white',
      colors =c(
        colorRampPalette(colors = c("#0F2C67","#BAC2D3"))(3),
        colorRampPalette(colors = c("#F0BCBC","#CD1818"))(3)
      ),name="OR"
    )+
    guides(x="axis_nested")+
    font("xy.text", size = 11,color = "black")+
    theme(axis.ticks = element_blank(),axis.title.x=element_blank(),
          axis.text.x=element_blank(),axis.title.y=element_blank())
  
  
  library(aplot)
  pdf("./plot/myeloid/GSE149614/marker_heat_dot.pdf", height = 7.5,width = 6.2)
  my_heatmap %>% insert_bottom(my_dot,height = 0.08)
  dev.off()
}

## 7.4 marker gene ano ------------------------------------------------------------------
markers1 <- readRDS("./outdata/GSE149614/myeloid/markers1.rds")
top10 <- markers1 %>% group_by(cluster) %>% top_n(10, avg_log2FC)


marker2 <- markers1[markers1$cluster=="11", ]
# marker2$gene
library(dplyr)
fine1 <- as.character(Myeloid3$seurat_clusters)
fine<-recode(fine1,
             "0" = "M0_Mø_SEPP1",
             "1" = "M1_Mø_SPP1",
             "2" = "M2_Mø_C1QC",
             "3" = "M3_Mon_S100A8",
             "4" = "M4_Mø_C1QB",
             "5" = "M5_Mon_FCN1",
             "6" = "M6_Mø_TOP2A",
             "7" = "M7_DC_CD1C",
             "8" = "M8_Mø_CD3D",
             "9" = "M9_DC_CLEC9A",
             "10" = "M10_Mast_TPSAB1",
             "11" = "M11_DC_GZMB"
)
main<-recode(fine1,
             "0" = "Macrophage",
             "1" = "Macrophage",
             "2" = "Macrophage",
             "3" = "Monocyte",
             "4" = "Macrophage",
             "5" = "Monocyte",
             "6" = "Macrophage",
             "7" = "DC",
             "8" = "Macrophage",
             "9" = "DC",
             "10" = "Mast",
             "11" = "DC"
)
Myeloid4 <- Myeloid3
Myeloid4$fine<-fine
Myeloid4$main<-main
dittoDimPlot(Myeloid4, reduction.use = "umap", var = "fine",
             do.label=F,size =0.2,opacity =0.8,legend.size=8)
dittoDimPlot(Myeloid4, reduction.use = "umap", var = "main",
             do.label=F,size =0.2,opacity =0.8,legend.size=8)
saveRDS(Myeloid4, file = "./outdata/GSE149614/myeloid/Myeloid4.rds")

pdf("./plot/myeloid/GSE149614/my_ano_plot.pdf", height = 4,width = 6)
dittoDimPlot(Myeloid4, reduction.use = "umap", var = "fine",
             do.label=F,size =0.2,opacity =0.8,legend.size=8)
dev.off()


