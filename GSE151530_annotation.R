
library(Seurat)
library(data.table)
# 1. initial processing ---------------------------------------------------

GSE151530 <- Read10X(data.dir ="./initial_data/GSE151530/")
cli <- fread("./initial_data/GSE151530/HCC.txt",data.table = F)
rownames(cli) <- cli$Cell
GSE151530 <- GSE151530[,match(cli$Cell, colnames(GSE151530))]
GSE151530_se <- CreateSeuratObject(GSE151530, meta.data = cli)
saveRDS(GSE151530_se, file = "./outdata/GSE151530/GSE151530_se.rds")

# 2.The cell categories in the articles were extracted --------------------------------------------------------------------
Tumor <- subset(GSE151530_se, Type=="Malignant cells")
Tcells <- subset(GSE151530_se, Type=="T cells")
Tcells <- subset(Tcells, Sample %in% unique(Tumor$Sample))

B <- subset(GSE151530_se, Type=="B cells")
B <- subset(B, Sample %in% unique(Tumor$Sample))

TAMs <- subset(GSE151530_se, Type=="TAMs")
TAMs <- subset(TAMs, Sample %in% unique(Tumor$Sample))


saveRDS(Tumor, file = "./outdata/GSE151530/Tumor.rds")
saveRDS(Tcells, file = "./outdata/GSE151530/Tcells.rds")
saveRDS(TAMs, file = "./outdata/GSE151530/TAMs.rds")
saveRDS(B, file = "./outdata/GSE151530/B.rds")



# 3. sample subtype -------------------------------------------------------

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
  Tumor <- readRDS("./outdata/GSE151530/Tumor.rds")
  Tumor <- NormalizeData(Tumor)
  data<-Tumor_pair(Tumor,result,result3)
  data<-na.omit(data)
  subtype <- Tumor_class(Tumor,data,group.by = "Sample")
  saveRDS(subtype, file = "./pair/GSE151530_subtype.rds")
}

# 4. T+B clustering---------------------------------------------------------------------
## 4.1 T+B first QC---------------------------------------------------------------------
Tcell <- readRDS("./outdata/GSE151530/Tcells.rds")
B <- readRDS("E:/project/HCC_immune_scRNA/outdata/GSE151530/B.rds")
##合并
Tcell<-merge(Tcell,B)
Tcell <- CreateSeuratObject(Tcell@assays$RNA@counts, meta.data = Tcell@meta.data)
Tcell1 <- NormalizeData(Tcell) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
  #ScaleData(vars.to.regress=c("percent.mt","percent.ribo"))
  ScaleData() %>% 
  RunPCA(verbose = F, npcs = 50)

Tcell1 <- Tcell1 %>% RunHarmony("Sample", plot_convergence = F,
                                          max.iter.harmony=10)

Tcell2 <- Tcell1 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.7) %>% 
  identity()

VlnPlot(object = Tcell2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,
        group.by="RNA_snn_res.0.7")
DimPlot(Tcell2,reduction ="umap",label = T)
DimPlot(Tcell2,reduction ="umap",label = T,group.by = "Type")
FeaturePlot(Tcell2,reduction = "umap",features = c("KLRF1","CD4", "CD79A","CD8A"))
FeaturePlot(Tcell2,reduction = "umap",features = c("CD3D","HBA1","IRF7"))
FeaturePlot(Tcell2,reduction = "umap",features = c("CD4","NKG7","CD3D"))
dittoDimPlot(Tcell2, reduction.use = "umap", var = "RNA_snn_res.0.5",
             do.label=T,size =0.2,opacity =0.8,legend.size=8)
#12低质量，11红细胞，13pDC
##提取
pDC <- subset(Tcell2, RNA_snn_res.0.7 %in% c("13"))
saveRDS(pDC, file = "./outdata/GSE151530/pDC.rds")
Ly <- subset(Tcell2, RNA_snn_res.0.7 %in% setdiff(as.character(0:13),c("11","12","13")))
saveRDS(Ly, file = "./outdata/GSE151530/Ly.rds")


## 4.2 T+B clustering -----------------------------------------------------------------
Ly <- readRDS("E:/project/HCC_immune_scRNA/outdata/GSE151530/Ly.rds")
Ly <- CreateSeuratObject(Ly@assays$RNA@counts, meta.data = Ly@meta.data)
Ly1 <- NormalizeData(Ly) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
  #ScaleData(vars.to.regress=c("percent.mt","percent.ribo"))
  ScaleData() %>% 
  RunPCA(verbose = F, npcs = 50)

Ly1 <- Ly1 %>% RunHarmony("Sample", plot_convergence = F,
                                max.iter.harmony=10)

Ly2 <- Ly1 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DimPlot(Ly2,reduction ="umap",label = T)
FeaturePlot(Ly2,reduction = "umap",features = c("CD79A","CD4","CD8A","KLRF1"))
VlnPlot(object = Ly2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,
        group.by="RNA_snn_res.0.5")

lymarker <- FindAllMarkers(Ly2, only.pos = T)
lymarker1 <- lymarker[lymarker$cluster=="0",]

##剔除群4，低质量且为双细胞
Ly3 <- subset(Ly2, RNA_snn_res.0.5 %in% c("4"))
DimPlot(Ly3,reduction ="umap",label = T)
FeaturePlot(Ly3,reduction = "umap",features = c("APOC1"))

Ly3 <- subset(Ly2, RNA_snn_res.0.5 %in% setdiff(as.character(0:9),c("4")))
lymarker <- FindAllMarkers(Ly3, only.pos = T)
lymarker1 <- lymarker[lymarker$cluster=="0",]
DimPlot(Ly3,reduction ="umap",label = T)

## 4.3 marker gene plot ------------------------------------------------------------------
### 4.3.1 heatmap --------------------------------------------------------------------

#平均热图(论文出图版)
{
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
  cluster_type <- data.frame(cluster = c("L0","L5","L7","L1","L2","L8","L3","L6","L9"), 
                             cell_type=c(rep("CD4+T", 3), rep("CD8+T", 3), rep("NK", 1), rep("B", 2)) )
  
  ##mean scale
  lymean <- AverageExpression(Ly3,features = marker_genes, group.by="seurat_clusters")
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
    geom_vline(xintercept = c(3.5, 6.5, 7.5), size = 1.5, color = "white") +
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

### 4.3.2 OR plot --------------------------------------------------------------------

##分类情况
subtype <- readRDS("./pair/GSE151530_subtype.rds") 

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
  
  metadata1 <- Ly3@meta.data
  temp <- metadata1[,c("Sample","seurat_clusters"),drop=F]
  temp$seurat_clusters <- paste0("L",temp$seurat_clusters)
  temp$subtype <- subtype[match(temp$Sample, rownames(subtype)), 4]
  
  library(data.table)
  library(plyr)
  cell_type<-data.table(temp[,-1])
  a1<-cell_type$seurat_clusters
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
  cluster_type <- data.frame(cluster = c("L0","L5","L7","L1","L2","L8","L3","L6","L9"), 
                             cell_type=c(rep("CD4+T", 3), rep("CD8+T", 3), rep("NK", 1), rep("B", 2)))
  ly1$celltype <- cluster_type$cell_type[match(ly1$rid, cluster_type$cluster)]
  ly1$celltype <- factor(ly1$celltype, levels = c("CD4+T", "CD8+T", "NK", "B"))
  
  ly_dot <- ggplot(ly1,aes(x = interaction(rid,celltype),y = cid,color = dirOR,size=plabel)) +
    geom_point() +
    scale_size_discrete(range=c(1,6))+
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
  pdf("./plot/lymphocyte/GSE151530/marker_heat_dot.pdf", height = 8,width = 6.3)
  ly_heatmap %>% insert_bottom(ly_dot,height = 0.08)
  dev.off()
}


####注释
main <- as.character(Ly3$RNA_snn_res.0.5)
main[main %in% c("5")]="Treg"
main[main %in% c("3")]="NK"
main[main %in% c("6","9")]="B"
main[main %in% c("0","7")]="CD4"
main[main %in% c("1","2","8")]="CD8"
Ly3$main <- main
DimPlot(Ly3,reduction ="umap",label = T,group.by = "main")
saveRDS(Ly3, file = "./outdata/GSE151530/Ly3.rds")

# 5. Myeloid clustering -----------------------------------------------------------------

## 5.1 QC ----------------------------------------------------------------
library(harmony)
TAMs <- readRDS("./outdata/GSE151530/TAMs.rds")
pDC <- readRDS("E:/project/HCC_immune_scRNA/outdata/GSE151530/pDC.rds")
##合并
TAM<-merge(TAMs,pDC)
TAM <- CreateSeuratObject(TAM@assays$RNA@counts, meta.data = TAM@meta.data)
TAM1 <- NormalizeData(TAM) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
  #ScaleData(vars.to.regress=c("percent.mt","percent.ribo"))
  ScaleData() %>% 
  RunPCA(verbose = F, npcs = 50)

TAM1 <- TAM1 %>% RunHarmony("Sample", plot_convergence = F,
                            max.iter.harmony=10)

TAM2 <- TAM1 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
DimPlot(TAM2,reduction ="umap",label = T)
DimPlot(TAM2,reduction ="umap",group.by = "Sample")

TAM2[["percent.mt"]] <- PercentageFeatureSet(TAM2, pattern = "^MT-")
TAM2[["nCount_RNA"]] = colSums(x = TAM2, slot = "counts")  # nCount_RNA
TAM2[["nFeature_RNA"]] = colSums(x = GetAssayData(object = TAM2, slot = "counts") > 0)
VlnPlot(object = TAM2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,
        group.by="RNA_snn_res.0.5")
mymarker0 <- FindAllMarkers(TAM2, only.pos = T)
mymarker<-mymarker0

AverageHeatmap(
  object = TAM2,                          # 分析对象
  markerGene = marker_genes,                 # 标记基因
  group.by = "RNA_snn_res.0.5",                # 按照聚类结果进行分组
  #cluster.order = as.character(c(0,2,3,4,11,13, 1,5,8,7,9,6,10,12)), # 聚类顺序
  cluster_columns = T,                   # 列不进行聚类
  #column_split = c(rep("CD4+T", 6), rep("CD8+T", 3), rep("NK", 2), rep("B", 3)), # 列分组
  row_split = c(rep("a", 4), rep("b", 5), rep("c", 5), rep("d", 3), rep("e", 4), rep("f", 8),
                rep("g", 4), rep("h", 2), rep("i", 3), rep("j", 3)), # 行分组
  border = TRUE                              # 边框
)
# 过滤核糖体
mymarker <- mymarker[!grepl("^RP[ SL ]", mymarker$gene, ignore.case = F),]
# 过滤线粒体
mymarker <- mymarker[!grepl("^MT-", mymarker$gene, ignore.case = F),]
mymarker1 <- mymarker[mymarker$cluster=="1",]

##剔除1群低质量且混杂细胞
Myeloid <- subset(TAM2, RNA_snn_res.0.5 %in% c("0",as.character(2:14)))
saveRDS(Myeloid, file = "./outdata/GSE151530/Myeloid.rds")

## 5.2 Myeloid reclustering -------------------------------------------------------------------
Myeloid<- readRDS("./outdata/GSE151530/Myeloid.rds")
Myeloid<-Myeloid %>% RunUMAP(reduction = "harmony", dims = 1:20,spread =0.8)
DimPlot(Myeloid,reduction ="umap",label = T)
FeaturePlot(Myeloid,reduction = "umap",features = c("C1QA"))
saveRDS(Myeloid,file = "./outdata/GSE151530/Myeloid2.rds")

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
# 生成矩阵气泡图，展示各类型巨噬细胞子集的表达情况
DotPlot(Myeloid, features = mac_name,
        assay='RNA' ,group.by = 'seurat_clusters') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
VlnPlot(object = Myeloid3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
marker<-FindAllMarkers(Myeloid,only.pos = T)
a<-marker[marker$cluster=="3",]

## 5.3 annotation plot ----------------------------------------------------------------
### 5.3.1 heatmap --------------------------------------------------------------
Myeloid <- readRDS("./outdata/GSE151530/Myeloid2.rds")
{
  cluster_type <- data.frame(cluster =  paste0("M", as.character(c(0,3,6,9,7,11,12,2,4,8,5,10,13,14))), cell_type=c(rep("Mø", 7), rep("Mon", 3), rep("DC", 3), rep("Mast", 1)))
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
  plotdata$variable <- factor(plotdata$variable, levels = paste0("M", as.character(c(0,3,6,9,7,11,12,2,4,8,5,10,13,14))))
  
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
    geom_vline(xintercept = c(7.5, 10.5, 13.5), size = 1.5, color = "white") +
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


### 5.3.2 OR ratio plot --------------------------------------------------------------
##分类情况
subtype <- readRDS("./pair/GSE151530_subtype.rds") 

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
  
  metadata1 <- Myeloid@meta.data
  temp <- metadata1[,c("Sample","seurat_clusters"),drop=F]
  temp$seurat_clusters <- paste0("M",temp$seurat_clusters)
  temp$subtype <- subtype[match(temp$Sample, rownames(subtype)), 4]
  
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
  cluster_type <- data.frame(cluster =  paste0("M", as.character(c(0,3,6,9,7,11,12,2,4,8,5,10,13,14))), cell_type=c(rep("Mø", 7), rep("Mon", 3), rep("DC", 3), rep("Mast", 1)))
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
  pdf("./plot/myeloid/GSE151530/marker_heat_dot.pdf", height = 7.5,width = 7)
  my_heatmap %>% insert_bottom(my_dot,height = 0.08)
  dev.off()
}

## 5.4 Myeloid annotation --------------------------------------------------------------
Myeloid2<-Myeloid
main <- as.character(Myeloid2$RNA_snn_res.0.5)
main[main %in% c("0","3","6","9","7","11","12")]="Macrophage"
main[main %in% c("4","2","8")]="Monocyte"
main[main %in% c("14")]="Mast"
main[main %in% c("5","10","13")]="DC"
Myeloid2$main <- main
DimPlot(Myeloid2,reduction ="umap",label = T,group.by = "main")
saveRDS(Myeloid2, file = "./outdata/GSE151530/Myeloid2.rds")
