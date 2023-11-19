

# 1. Construct a single cell model -----------------------------------------------------------------
## 1.1 Four gene pair groups ---------------------------------------------------------------------
## 2000 HVFGs from GSE149614(Cancer cells)
library(Seurat)
Tumor <- readRDS("./outdata/GSE149614/Tumor.rds")
Tumor <- NormalizeData(Tumor) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

##HVGs in 4 groups
vst<-Tumor@assays[["RNA"]]@var.features
gene_group <- list()
gene_group$group1 <- vst[1:500]
gene_group$group2 <- vst[501:1000]
gene_group$group3 <- vst[1001:1500]
gene_group$group4 <- vst[1501:2000]
saveRDS(gene_group, file = "./pair/gene_group.rds")

## bulk data gene pairs
Bulk_exp <- readRDS("./data/Bulk_exp.rds")
pair1 <- lapply(gene_group, function(x){
  gene_name <- intersect(x,Bulk_exp$GSE14520$Symbol)
  gene_name <- intersect(gene_name,Bulk_exp$TCGA$SYMBOL)
  gene_name <- intersect(gene_name,Bulk_exp$ICGC$symbol)
  gene_name <- intersect(gene_name,Bulk_exp$GSE116174$Symbol)
  gene_name <- intersect(gene_name,Bulk_exp$`In-house`$V1)
  pair_combn <-t(combn(gene_name, 2))
  return(pair_combn)
})
saveRDS(pair1, file = "./pair/pair1.rds")

##得到每个分组的稳定基因对,p<0.05,稳定对阈值=0.7
library(purrr)
predict_data <- readRDS("./data/predict_data.rds")
pair_fisher <- function(x,y, pair, subtype){
  data <- egg_us(x)
  data <- data[pair[,1],] - data[pair[,2],]
  data[]<-apply(data,2,function(x) ifelse(x>0,1,x))
  data[]<-apply(data,2,function(x) ifelse(x<0,-1,x))
  #
  print(111)
  out1<-apply(data,1,function(x){
    a=sum(y$Subtype[which(x==1)]==subtype[1])
    b=sum(y$Subtype[which(x==1)]==subtype[2])
    d=sum(y$Subtype[which(x==-1)]==subtype[1])
    e=sum(y$Subtype[which(x==-1)]==subtype[2])
    cbind(a,b,d,e)
  })
  out1<-t(out1)
  print(222)
  pvalue<-apply(out1,1,function(x){
    estimate<-fisher.test(matrix(c(x[1],x[2],x[3],x[4]),ncol=2))$estimate
    if(estimate<1){
      out<-cbind(fisher.test(matrix(c(x[3],x[4],x[1],x[2]),ncol=2))$p.value,-1,
                 fisher.test(matrix(c(x[3],x[4],x[1],x[2]),ncol=2))$estimate,
                 x[3],x[4],x[1],x[2],x[3]/(x[1] + x[3]),  x[2]/(x[2] + x[4]))
    }else{
      out<-cbind(fisher.test(matrix(c(x[1],x[2],x[3],x[4]),ncol=2))$p.value,1,
                 fisher.test(matrix(c(x[1],x[2],x[3],x[4]),ncol=2))$estimate,
                 x[1],x[2],x[3],x[4],x[1]/(x[1] + x[3]), x[4]/(x[2] + x[4]))
    }
  })
  pvalue<-t(pvalue)
  out<-data.frame(pair,pvalue)
  out[which(out[,4]<0),1:2] <- out[which(out[,4]<0),2:1]
  colnames(out)[3]<-"p"
  out$name <- paste0(out[,1],"_",out[,2])
  return(out)
}
egg_us<-function(data){
  rownames(data)<-data[,1]
  data<-data[,-1]
  data<-as.matrix(data)
  return(data)
}
result <-lapply(pair1, function(o){
  print("<<<1>>>>")
  pair<-o
  bulk_pair <- list()
  bulk_pair1 <- list()
  #TCGA
  bulk_pair$TCGA1 <- pair_fisher(Bulk_exp$TCGA,predict_data$TCGA,pair[,-3],c("C1","C2"))
  bulk_pair1$TCGA2 <- pair_fisher(Bulk_exp$TCGA,predict_data$TCGA,pair[,-3],c("C2","C1"))
  #ICGC
  bulk_pair$ICGC1 <- pair_fisher(Bulk_exp$ICGC,predict_data$ICGC,pair[,-3],c("C1","C2"))
  bulk_pair1$ICGC2 <- pair_fisher(Bulk_exp$ICGC,predict_data$ICGC,pair[,-3],c("C2","C1"))
  
  #GSE14520
  bulk_pair$GSE145201 <- pair_fisher(Bulk_exp$GSE14520,predict_data$GSE14520,pair[,-3],c("C1","C2"))
  bulk_pair1$GSE145202 <- pair_fisher(Bulk_exp$GSE14520,predict_data$GSE14520,pair[,-3],c("C2","C1"))
  
  #GSE116174
  bulk_pair$GSE1161741 <- pair_fisher(Bulk_exp$GSE116174,predict_data$GSE116174,pair[,-3],c("C1","C2"))
  bulk_pair1$GSE1161742 <- pair_fisher(Bulk_exp$GSE116174,predict_data$GSE116174,pair[,-3],c("C2","C1"))
  
  #Inhouse
  bulk_pair$Inhouse1 <- pair_fisher(Bulk_exp$`In-house`,predict_data$`In-house`,pair[,-3],c("C1","C2"))
  bulk_pair1$Inhouse2 <- pair_fisher(Bulk_exp$`In-house`,predict_data$`In-house`,pair[,-3],c("C2","C1"))
  
  pair_bulk_fil <- lapply(bulk_pair, function(x){
    test <- x[x$p < 0.05 & x$X8>0.7 ,] 
    return(test$name)
  })
  pair_temp <- reduce(pair_bulk_fil,intersect)
  pair_temp <- bulk_pair$TCGA[match(pair_temp, bulk_pair$TCGA$name),1:2]
  
  ###
  pair_bulk_fil1 <- lapply(bulk_pair1, function(x){
    test <- x[x$p < 0.05 & x$X8>0.7 ,] 
    return(test$name)
  })
  pair_temp1 <- reduce(pair_bulk_fil1,intersect)
  pair_temp1 <- bulk_pair1$TCGA[match(pair_temp1, bulk_pair1$TCGA$name),1:2]
  
  
  test2 <- mapply(function(k,j){
    data <- egg_us(k)
    data <- data[pair_temp[,1],] - data[pair_temp[,2],]
    out<-apply(data, 2, function(x){
      sum(x>0)/length(x)
    })
    
    data <- egg_us(k)
    data <- data[pair_temp1[,1],] - data[pair_temp1[,2],]
    out1<-apply(data, 2, function(x){
      sum(x>0)/length(x)
    })
    
    test1<-data.frame(C1=out,C2=out1,check.names = F)
  },Bulk_exp,SIMPLIFY = F)
  print("<<<2>>>>")
  test2$pair_temp <- pair_temp
  test2$pair_temp1 <- pair_temp1
  return(test2)
})

saveRDS(result, file = "./pair/result.rds")

result3 <- list()
name="TCGA"
result3[[name]] <- lapply(result, function(x){
  x[[name]]
})
result3[[name]] <- do.call(cbind,lapply(result3[[name]], function(x) x))
result3[[name]] <- data.frame(subtype=predict_data[[name]]$Subtype,result3[[name]])

name="ICGC"
result3[[name]] <- lapply(result, function(x){
  x[[name]]
})
result3[[name]] <- do.call(cbind,lapply(result3[[name]], function(x) x))
result3[[name]] <- data.frame(subtype=predict_data[[name]]$Subtype,result3[[name]])

name="GSE14520"
result3[[name]] <- lapply(result, function(x){
  x[[name]]
})
result3[[name]] <- do.call(cbind,lapply(result3[[name]], function(x) x))
result3[[name]] <- data.frame(subtype=predict_data[[name]]$Subtype,result3[[name]])

name="GSE116174"
result3[[name]] <- lapply(result, function(x){
  x[[name]]
})
result3[[name]] <- do.call(cbind,lapply(result3[[name]], function(x) x))
result3[[name]] <- data.frame(subtype=predict_data[[name]]$Subtype,result3[[name]])

name="In-house"
result3[[name]] <- lapply(result, function(x){
  x[[name]]
})
result3[[name]] <- do.call(cbind,lapply(result3[[name]], function(x) x))
result3[[name]] <- data.frame(subtype=predict_data[[name]]$Subtype,result3[[name]])
saveRDS(result3, file = "./pair/result3.rds")

## 1.2 knn ---------------------------------------------------------------------
library(caret)
set.seed(456)
f = 10  # f folds resampling
r = 3 # r repeats
train<-rbind(result3$TCGA, result3$GSE14520)
ctrl <- trainControl(method="repeatedcv",
                     number = f, ## 10-folds cv
                     classProbs=TRUE,
                     repeats = r, ## 3-repeats cv,
)
kknn = expand.grid(kmax = c(5,7,9,11,13,15,17,19,21,25,31,37,45), distance = c(2) , kernel = 'optimal')
resultknn3<-train(subtype ~ .,
                  data = train,
                  method = "kknn",
                  metric="Accuracy",
                  trControl=ctrl,
                  tuneGrid = kknn)
saveRDS(resultknn3, file = "./pair/resultknn3.rds")

##bulk数据集验证
  name=c("TCGA", "ICGC", "GSE14520", "GSE116174", "In-house")
  acc2<-lapply(name, function(x){
    a <- predict(resultknn3, newdata = result3[[x]])
    aa<-result3[[x]]$subtype
    aa<-factor(aa,levels=c("C1","C2"))
    #print(confusionMatrix(a,aa)[["byClass"]][["Balanced Accuracy"]])
    #confusionMatrix(a,aa)[["byClass"]][["F1"]]
    confusionMatrix(a,aa)[["overall"]][["Accuracy"]]
  })
  acc2<-unlist(acc2)
  acc2
  
## 1.3 nb ----------------------------------------------------------------------
  set.seed(123)
  nb = expand.grid(fL =  c(0,0.5,1,1.5,2.0), usekernel = TRUE, adjust = c(0.5,0.75,1,1.25,1.5))  
  nb <- train(subtype ~ .,
              data = train,
              method = "nb",
              metric="Accuracy",
              trControl=ctrl,
              tuneGrid = nb)
  saveRDS(nb, file = "./pair/nb.rds")
  
##bulk数据集验证
  name=c("TCGA", "ICGC", "GSE14520", "GSE116174", "In-house")
  acc1<-lapply(name, function(x){
    a <- predict(nb, newdata = result3[[x]])
    aa<-result3[[x]]$subtype
    aa<-factor(aa,levels=c("C1","C2"))
    #print(confusionMatrix(a,aa)[["byClass"]][["Balanced Accuracy"]])
    #confusionMatrix(a,aa)[["byClass"]][["F1"]]
    confusionMatrix(a,aa)[["overall"]][["Accuracy"]]
  })
  acc1<-unlist(acc1)
  acc1

## 1.4 Assessment of single-cell classification models in bulk transcriptomic datasets---------------------------------------------------------------------

### 1.4.1 point cor -------------------------------------------------------------------
  #knn
  {
  KNN <- readRDS("./pair/0.7/resultknn3.rds")
  result3 <- readRDS("./pair/result3.rds")
  predict_data <- readRDS("./data/predict_data.rds")
  ##  prob
  name=c("TCGA", "ICGC", "GSE14520", "GSE116174", "In-house")
  prob<-lapply(name, function(x){
    a <- predict(KNN, newdata = result3[[x]],type="prob")
    a1 <- data.frame(a,model=predict(KNN, newdata = result3[[x]]),
                     model_com = predict_data[[x]]$Subtype,
                     Subtype = predict_data[[x]]$Subtype,
                     Ratio = predict_data[[x]]$Ratio)
    a1$model_com[a1$model_com != a1$model] <- "Missubtype"
    a1
  })
  names(prob)<-name
  
  ##plot
  ##C1 and C2
  library(dplyr)
  library(ggplot2)
  library(ggstatsplot)
  library(rlang)
  library(tidyr)
  library(ggpubr)
  ##The scatter plot of distribution of the two scaled-scores for all samples in four datasets.
  KNN_plot <- lapply(name, function(x){
    print(x)
    Consistency <- (table(prob[[x]]$model_com)[1]+table(prob[[x]]$model_com)[2])/length(prob[[x]]$model_com)
    cor_result<-cor.test(prob[[x]]$Ratio,prob[[x]]$C1,method="spearman")
    if(cor_result$p.value < 0.001){
      cortext<-paste0("R = ",round(cor_result$estimate,2)," ,p < 0.001")
    }
    p1<-point_section(prob[[x]],x=Ratio,y=C1,
                      colour=model_com,colour_value=c("#CD1818", "#0F2C67", "grey70"),
                      annotate="",
                      xprobs = 0.5,
                      intercept = 0.5,
                      slope=0,
                      xintercept=0.5,
                      yintercept=0.5,
                      cor=F,
                      sizeP = 2,
                      x_label = "The risk for C1 by 29 pairs",
                      y_label = "KNN probability",
                      col_section_up = "#4DBBD5FF",
                      col_section_down = "#E64B35FF")+ ggtitle(x)
    p1+
      annotate('text',
               x=0,
               y=0.89,
               label= paste0("Consistency = ",round(Consistency,2)),
               size=5,color='black',hjust = "left")+
      annotate('text',
               x=0,
               y=0.99,
               label= cortext,
               size=5,color='black',hjust = "left")
  })
  names(KNN_plot) <- name
}

  #nb
  {
    nb <- readRDS("./pair/nb.rds")
    result3 <- readRDS("./pair/result3.rds")
    predict_data <- readRDS("./data/predict_data.rds")
    ##  prob
    name=c("TCGA", "ICGC", "GSE14520", "GSE116174", "In-house")
    prob<-lapply(name, function(x){
      a <- predict(nb, newdata = result3[[x]],type="prob")
      a1 <- data.frame(a,model=predict(nb, newdata = result3[[x]]),
                       model_com = predict_data[[x]]$Subtype,
                       Subtype = predict_data[[x]]$Subtype,
                       Ratio = predict_data[[x]]$Ratio)
      a1$model_com[a1$model_com != a1$model] <- "Missubtype"
      a1
    })
    names(prob)<-name
    
    ##plot
    ##C1 and C2
    library(dplyr)
    library(ggplot2)
    library(ggstatsplot)
    library(rlang)
    library(tidyr)
    library(ggpubr)
    ##The scatter plot of distribution of the two scaled-scores for all samples in four datasets.
    nb_plot <- lapply(name, function(x){
      print(x)
      Consistency <- (table(prob[[x]]$model_com)[1]+table(prob[[x]]$model_com)[2])/length(prob[[x]]$model_com)
      cor_result<-cor.test(prob[[x]]$Ratio,prob[[x]]$C1,method="spearman")
      if(cor_result$p.value < 0.001){
        cortext<-paste0("R = ",round(cor_result$estimate,2)," ,p < 0.001")
      }
      p1<-point_section(prob[[x]],x=Ratio,y=C1,
                        colour=model_com,colour_value=c("#CD1818", "#0F2C67", "grey70"),
                        annotate="",
                        xprobs = 0.5,
                        intercept = 0.5,
                        slope=0,
                        xintercept=0.5,
                        yintercept=0.5,
                        cor=F,
                        sizeP = 2,
                        x_label = "The risk for C1 by 29 pairs",
                        y_label = "NB probability",
                        col_section_up = "#4DBBD5FF",
                        col_section_down = "#E64B35FF")+ ggtitle(x)
      p1+
        annotate('text',
                 x=0,
                 y=0.89,
                 label= paste0("Consistency = ",round(Consistency,2)),
                 size=5,color='black',hjust = "left")+
        annotate('text',
                 x=0,
                 y=0.99,
                 label= cortext,
                 size=5,color='black',hjust = "left")
    })
    names(nb_plot) <- name
  }  

#mean
  {
    result3 <- readRDS("./pair/result3.rds")
    predict_data <- readRDS("./data/predict_data.rds")
    ##  prob
    name=c("TCGA", "ICGC", "GSE14520", "GSE116174", "In-house")
    prob<-lapply(name, function(x){
      a <-data.frame(C1=apply(result3[[x]][,c(2,4,6,8)],1,mean),
                     C2=apply(result3[[x]][,c(3,5,7,9)],1,mean))
      temp <- ifelse(a$C1>a$C2,"C1","C2")
      a1 <- data.frame(a,model=temp,
                       model_com = predict_data[[x]]$Subtype,
                       Subtype = predict_data[[x]]$Subtype,
                       Ratio = predict_data[[x]]$Ratio,
                       meandif=a$C1-a$C2)
      a1$model_com[a1$model_com != a1$model] <- "Missubtype"
      a1
    })
    names(prob)<-name
    
    ##plot
    ##C1 and C2
    library(dplyr)
    library(ggplot2)
    library(ggstatsplot)
    library(rlang)
    library(tidyr)
    library(ggpubr)
    ##The scatter plot of distribution of the two scaled-scores for all samples in four datasets.
    mean_plot <- lapply(name, function(x){
      print(x)
      Consistency <- (table(prob[[x]]$model_com)[1]+table(prob[[x]]$model_com)[2])/length(prob[[x]]$model_com)
      cor_result<-cor.test(prob[[x]]$Ratio,prob[[x]]$C1,method="spearman")
      if(cor_result$p.value < 0.001){
        cortext<-paste0("R = ",round(cor_result$estimate,2)," ,p < 0.001")
      }
      p1<-point_section(prob[[x]],x=Ratio,y=meandif,
                        colour=model_com,colour_value=c("#CD1818", "#0F2C67", "grey70"),
                        annotate="",
                        xprobs = 0.5,
                        intercept = 0,
                        slope=0,
                        xintercept=0.5,
                        yintercept=0,
                        cor=F,
                        sizeP = 2,
                        x_label = "The risk for C1 by 29 pairs",
                        y_label = "Mean probability",
                        col_section_up = "#4DBBD5FF",
                        col_section_down = "#E64B35FF")+ ggtitle(x)
      p1+
        annotate('text',
                 x=0,
                 y=0.89,
                 label= paste0("Consistency = ",round(Consistency,2)),
                 size=5,color='black',hjust = "left")+
        annotate('text',
                 x=0,
                 y=0.99,
                 label= cortext,
                 size=5,color='black',hjust = "left")
    })
    names(mean_plot) <- name
  }
  
  ##plot
  library(patchwork)
  pdf("./plot/bulk_point.pdf",height = 11,width = 20)
  (wrap_plots(KNN_plot,ncol = 5)+plot_layout(guides = 'collect'))/
    (wrap_plots(nb_plot,ncol = 5)+plot_layout(guides = 'collect'))/
    (wrap_plots(mean_plot,ncol = 5)+plot_layout(guides = 'collect'))
  dev.off()
  
### 1.4.2 heatmap -------------------------------------------------------------------
  library(caret)
  name=c("TCGA", "ICGC", "GSE14520", "GSE116174", "In-house")
  knn_acc <- lapply(name, function(x){
    acc = sum(KNN_plot[[x]]$data$model_com != "Missubtype") / length(KNN_plot[[x]]$data$model_com)
    return(acc)
  })
  knn_acc <- unlist(knn_acc)
  
  nb_acc <- lapply(name, function(x){
    acc = sum(nb_plot[[x]]$data$model_com != "Missubtype") / length(nb_plot[[x]]$data$model_com)
    return(acc)
  })
  nb_acc <- unlist(nb_acc)
  
  mean_acc <- lapply(name, function(x){
    acc = sum(mean_plot[[x]]$data$model_com != "Missubtype") / length(mean_plot[[x]]$data$model_com)
    return(acc)
  })
  mean_acc <- unlist(mean_acc)
  
  eggtest<-data.frame(KNN=knn_acc,NB=nb_acc,MEAN=mean_acc)
  rownames(eggtest)<-c("TCGA","ICGC","GSE14520","GSE116174","In-house")
  
  ###plot
  ####
  plotdata <- t(eggtest)
  library(ggsci)
  library(pheatmap)
  library("viridis")  
  name = c("TCGA","ICGC","GSE14520","GSE116174","In-house")
  name<-as.data.frame(name)
  rownames(name)<-c("TCGA","ICGC","GSE14520","GSE116174","In-house")
  colnames(name)<-"Dataset"
  ann_colors = list(
    Dataset = c("TCGA"  = pal_nejm("default")(5)[1],
                "ICGC" = pal_nejm("default")(5)[2], 
                "GSE14520" = pal_nejm("default")(5)[3],
                "GSE116174" = pal_nejm("default")(5)[5],
                "In-house" = pal_nejm("default")(5)[4])
  )
  
  p1<-pheatmap(plotdata[1:3,],cluster_cols=F,cluster_rows=F,breaks = c(seq(0.5,1,0.01)),
               color =   viridis(51,alpha =1, option = "D"),
               na_col = "white",annotation_col = name,gaps_col=1:5,show_colnames=F,
               annotation_colors=ann_colors,
               display_numbers=T,number_color="black",
               fontsize_number=13,border_color="white") 
  pdf("./plot/threemodel_bulk_acc.pdf",height = 2,width = 5)
  print(p1)
  dev.off()

## 1.5 single cell heatmap---------------------------------------------------------------------

### 1.5.1 GSE149614-------------------------------------------------------------------
  {
    library(pheatmap)
    library(GSVA)
    library(Seurat)
    library(dplyr)
    library(patchwork)
    ####
    ##sample subtype
    {
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
      source("./code/validation.R")
      Tumor <- readRDS("./outdata/GSE149614/Tumor.rds")
      data<-Tumor_pair(Tumor,result,result3)
      data<-na.omit(data)
      subtype <- Tumor_class(Tumor,data,group.by = "Sample")
      #saveRDS(subtype, file = "./pair/GSE149614_subtype.rds")
    }
    
    ##plot
    library(ComplexHeatmap)
    library(ggsci)
    library(circlize)
    library(dittoSeq)
    ###Hep_plot
    nejm <- dittoColors()[1:10]
    #nejm <- pal_npg("nrc")(10)
    names(nejm) <- rownames(subtype)
    column_ha = HeatmapAnnotation(Sample=rownames(subtype),
                                  "KNN"=subtype$knn,
                                  "NB"= subtype$nb,
                                  "Mean" = subtype$aver,
                                  "Final Subtype" = subtype$subtype,
                                  col  = list(
                                    Sample=nejm,
                                    "KNN"= c("C1" = "#CD1818", "C2" = "#0F2C67"),
                                    "NB" = c("C1" = "#CD1818", "C2" = "#0F2C67"),
                                    "Mean" = c("C1" = "#CD1818", "C2" = "#0F2C67"),
                                    "Final Subtype" = c("C1" = "#CD1818", "C2" = "#0F2C67")
                                  ),annotation_name_side="left")
    temp <- matrix(0,nrow = 1,ncol = 10)
    
    p1<-Heatmap(temp,col="white",cluster_columns = F,cluster_rows = F,show_row_names=F,show_column_names=F,
                top_annotation = column_ha,show_heatmap_legend=F)
    p2<-draw(p1,annotation_legend_side="bottom")
    pdf("./plot/hep_single_model_GSE149614.pdf",width=8,height=3.5)
    print(p2)
    dev.off()
  }
### 1.5.2 GSE151530 -------------------------------------------------------------------
  {
    library(pheatmap)
    library(GSVA)
    library(Seurat)
    library(dplyr)
    library(patchwork)
    ####
    ##sample subtype
    {
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
      source("./code/validation.R")
      Tumor <- readRDS("./outdata/GSE151530/Tumor.rds")
      Tumor <- NormalizeData(Tumor)
      data<-Tumor_pair(Tumor,result,result3)
      data<-na.omit(data)
      subtype <- Tumor_class(Tumor,data,group.by = "Sample")
    }
    
    ##plot
    library(ComplexHeatmap)
    library(ggsci)
    library(circlize)
    library(dittoSeq)
    ###Hep_plot
    nejm <- dittoColors()[1:20]
    #nejm <- pal_npg("nrc")(10)
    names(nejm) <- rownames(subtype)
    column_ha = HeatmapAnnotation(Sample=rownames(subtype),
                                  "KNN"=subtype$knn,
                                  "NB"= subtype$nb,
                                  "Mean" = subtype$aver,
                                  "Final Subtype" = subtype$subtype,
                                  col  = list(
                                    Sample=nejm,
                                    "KNN"= c("C1" = "#CD1818", "C2" = "#0F2C67"),
                                    "NB" = c("C1" = "#CD1818", "C2" = "#0F2C67"),
                                    "Mean" = c("C1" = "#CD1818", "C2" = "#0F2C67"),
                                    "Final Subtype" = c("C1" = "#CD1818", "C2" = "#0F2C67")
                                  ),annotation_name_side="left")
    temp <- matrix(0,nrow = 1,ncol = 20)
    
    p1<-Heatmap(temp,col="white",cluster_columns = F,cluster_rows = F,show_row_names=F,show_column_names=F,
                top_annotation = column_ha,show_heatmap_legend=F)
    p2<-draw(p1,annotation_legend_side="bottom")
    pdf("./plot/hep_single_model_GSE151530.pdf",width=8,height=5)
    print(p2)
    dev.off()
  }
  
  
# 2. Single cell subtype classification and functional mapping --------------------------------------------------------------------
  
## 2.1 load data ----------------------------------------------------------------
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
## 2.2 GSE149614 Tumor---------------------------------------------------------------
### 2.2.1  Predicted sample subtype --------------------------------------------------------------
  ##sample subtype
  {
  Tumor <- readRDS("./outdata/GSE149614/Tumor.rds")
  data<-Tumor_pair(Tumor,result,result3)
  data<-na.omit(data)
  subtype <- Tumor_class(Tumor,data,group.by = "Sample")
  saveRDS(subtype, file = "./pair/GSE149614_subtype.rds")
  }
  

### 2.2.2  Five malignant scores --------------------------------------------------------------
  ##肿瘤细胞计算得分分类并画图
  Tumor <- tumor_score(Tumor,subtype = subtype, group.by="Sample")
  saveRDS(Tumor, file = "./outdata/GSE149614/Tumor_subtype.rds")
  name=c("ISG.RS", "Proliferation","CTA","Stem","Hypoxia","Glycloysis",
         "SPP1","ITGA5","ITGB1")
  Tumormeta <- Tumor@meta.data[,c("subtype","seurat_clusters","Proliferation","CTA","Stem","Hypoxia","Glycloysis")]
  Tumormeta$seurat_clusters <- paste0("HCC", Tumormeta$seurat_clusters)
  Tumormeta <- reshape2::melt(Tumormeta,id.vars=c("subtype","seurat_clusters"))
  Tumormeta$seurat_clusters <- factor(Tumormeta$seurat_clusters, levels = c("HCC1","HCC2","HCC6","HCC9","HCC10","HCC14"))
  rm(Tumor)
  
  ##画图
  library(dittoSeq)
  library(patchwork)
  name <- as.character(unique(Tumormeta$variable))
  
  ##cluster
  {
  plot_cluster <- list() 
  for (i in 1:5) {
    plot_cluster[[name[i]]]<-plot.violin(Tumormeta[Tumormeta$variable == name[i],],x=seurat_clusters,y=value,
              fill = seurat_clusters,
              palette.fill=dittoColors()[c(1,2,6,9,10,14)+1],
              palette.color = dittoColors()[c(1,2,6,9,10,14)+1],
              color = "white",
              fill.b = "white",
              color.b = seurat_clusters,
              #group=seurat_clusters,
              alpha.v = 0.7,
              notch.b=F,
              width.b=0.5,
              lwd.b = 0.8,
              add = c("boxplot","violinplot"),#,"shadow"
              x_label="",
              y_label = name[i],
              theme = theme_bw(base_size = 16))+
    theme(axis.text.x = element_blank(),   ## 删去所有刻度标签
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 1,unit = 'cm'))+
    theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5))+
    coord_cartesian(clip = 'off')#+
    #theme(legend.position = 'none')
  }
  }
  
  ##cluster && subtype
  {
  plot_cluster_subtype <- list() 
  for (i in 1:5) {
    plot_cluster_subtype[[name[i]]]<-plot.violin(Tumormeta[Tumormeta$variable == name[i],],x=seurat_clusters,y=value,
                                           fill = subtype,
                                           palette.fill=c("#CD1818","#0F2C67"),
                                           palette.color = c("#CD1818","#0F2C67"),
                                           color = "white",
                                           fill.b = "white",
                                           color.b = subtype,
                                           group=subtype,
                                           alpha.v = 0.7,
                                           notch.b=F,
                                           width.b=0.5,
                                           add = c("boxplot","violinplot","shadow"),#,"shadow"
                                           lwd.b = 0.8,
                                           shadow.col = c("grey40","white"),
                                           x_label="",
                                           y_label = "",
                                           theme = theme_bw(base_size = 16))+
      theme(axis.text.x = element_blank(),   ## 删去所有刻度标签
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),   ## 删去所有刻度标签
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 1,unit = 'cm'))+
      theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5))+
      coord_cartesian(clip = 'off')+
      theme(legend.position = 'none')
  }
  }
 
  ##subtype
  {
  plot_subtype <- list() 
  for (i in 1:5) {
    plot_subtype[[name[i]]]<-plot.violin(Tumormeta[Tumormeta$variable == name[i],],x=name[i],y=value,
                                           fill = subtype,
                                           palette.fill=c("#CD1818","#0F2C67"),
                                           palette.color = c("#CD1818","#0F2C67"),
                                           color = "white",
                                           fill.b = "white",
                                           color.b = subtype,
                                           group=subtype,
                                           alpha.v = 0.7,
                                           notch.b=F,
                                           width.b=0.5,
                                           add = c("boxplot","violinplot"),#,"shadow"
                                         lwd.b = 0.8,
                                           x_label="",
                                           y_label = "",
                                           theme = theme_bw(base_size = 16))+
      theme(axis.text.x = element_blank(),   ## 删去所有刻度标签
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),   ## 删去所有刻度标签
            axis.ticks.y = element_blank(),
            plot.margin = margin(t = 1,unit = 'cm'))+
      theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5))+
      coord_cartesian(clip = 'off')
  }
  }
  
  #plot
  {
  library(jjAnno)
  library(dittoSeq)
  plot_CTA <- annoRect(object = plot_cluster$CTA,
                       annoPos = 'top',
                       xPosition = c(1:6),
                       pFill = dittoColors()[c(1,2,6,9,10,14)+1],
                       pCol = dittoColors()[c(1,2,6,9,10,14)+1],
                       yPosition = c(-0.33,-0.38),
                       rectWidth = 1)
  plot_CTA_cs <- annoRect(object = plot_cluster_subtype$CTA,
                       annoPos = 'top',
                       xPosition = c(1:6),
                       pFill = dittoColors()[c(1,2,6,9,10,14)+1],
                       pCol = dittoColors()[c(1,2,6,9,10,14)+1],
                       yPosition = c(-0.33,-0.38),
                       rectWidth = 1)
  
  pdf(file = "./plot/cancer/malignant.pdf",width =10 ,height = 9)
  (plot_cluster$Glycloysis + plot_cluster_subtype$Glycloysis + plot_subtype$Glycloysis+plot_layout(widths = c(2,4,1))) / 
  (plot_cluster$Hypoxia + plot_cluster_subtype$Hypoxia + plot_subtype$Hypoxia+plot_layout(widths = c(2,4,1))) / 
  (plot_cluster$Proliferation + plot_cluster_subtype$Proliferation + plot_subtype$Proliferation+plot_layout(widths = c(2,4,1))) / 
  (plot_cluster$Stem + plot_cluster_subtype$Stem + plot_subtype$Stem+plot_layout(widths = c(2,4,1))) / 
  (plot_CTA + plot_CTA_cs + plot_subtype$CTA+plot_layout(widths = c(2,4,1))) + plot_layout(guides = "collect")
  dev.off()
  }

  
### 2.2.3  barplot & OR  --------------------------------------------------------------  
  ##barplot
  {
  library(dittoSeq)
  Tumor <- readRDS("./outdata/GSE149614/Tumor_subtype.rds")
  Tumor$clusters <- paste0("HCC",Tumor$seurat_clusters)
  Tumor$clusters <- factor(Tumor$clusters, levels = c("HCC1","HCC2","HCC6","HCC9","HCC10","HCC14"))
  p1 <- dittoBarPlot(Tumor, var = "clusters", group.by = "subtype", main = "Cancer",retain.factor.levels=T,
                     color.panel=dittoColors()[c(1,2,6,9,10,14)+1])+
    theme(axis.title.x =element_text(size=13,color = "black"), 
          axis.title.y=element_text(size=13,color = "black"),
          axis.text.x=element_text(size=13,color = "black"),
          axis.text.y=element_text(size=13,color = "black"))
  }
  
  ##pie plot
  {
  p2 <- dittoBarPlot(Tumor, var = "subtype", group.by = "clusters", main = "Cancer")$data
  out1 <- p2[p2$grouping=="HCC14", ]
  out2 <- out1 %>% mutate(lab.pos = cumsum(percent*100) - 0.5*100*percent)
  library(ggsci)
  library(ggh4x)

  p2 <- ggplot(out2, aes(x = "", y = percent*100, fill = label)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0)+
    geom_text(aes(y = lab.pos, label = round(percent*100, 2)), color = "white",size=5)+
    theme_void()+
    scale_fill_manual(values=c("#CD1818","#0F2C67"))+
    scale_colour_manual(values=c("#CD1818","#0F2C67"))
  }
 
  
  ##OR
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
  
  {
  metadata1 <- Tumor@meta.data
  temp <- metadata1[,c("subtype","clusters")]
  temp <- temp[order(temp$clusters),]
  library(data.table)
  library(plyr)
  cell_type<-data.table(temp)
  a1<-cell_type$clusters
  a2<-cell_type$subtype
  temp <- unclass(cell_type[,table(a1,a2)])
  temp <- temp[unique(a1),]
  my<-test.dist.table(temp)
  my$OR <- ifelse(my$OR>4, 4,my$OR)
  my$rid <- factor(my$rid, levels = c("HCC1","HCC2","HCC6","HCC9","HCC10","HCC14"))
  
  p3 <- ggplot(my,aes(x=cid,y=rid,fill=OR))+
    geom_raster()+
    scale_fill_viridis_c( option = "D",alpha =1,begin = 0, end = 1)+
    xlab(NULL) + ylab(NULL)+
    theme(panel.grid = element_blank(),panel.background =element_blank() )+
    geom_tile()+
    theme(axis.ticks.x = element_blank(),axis.ticks.y=element_blank())+
    font("xy", size = 11)+
    font("xy.text", size = 11,color = "black")+
    font("legend.text", size = 12)+
    font("legend.title", size = 12)+
    theme(axis.title.x =element_blank())
  }

  
  pdf(file = "./plot/cancer/percent_OR.pdf",width =6.5 ,height = 2.5)
  p1 + p2 + p3 + plot_layout(widths =c(1,2.5,1))
  dev.off()
  
  
### 2.2.4  Marker genes ( abundance, malignant, survival) --------------------------------------------------------------
  library(GSVA)
  library(ggpubr)
  library(ggh4x)
  Bulk_exp <- readRDS("./data/Bulk_exp.rds")
  predict_data <- readRDS("./data/predict_data.rds")
  ##marker gene
  marker<-FindAllMarkers(Tumor,only.pos = T,min.pct=0.5)
  markers1 <- marker[!grepl("^RP[ SL ]", marker$gene, ignore.case = F),]
  markers1 <- markers1[!grepl("^MT-", markers1$gene, ignore.case = F),]
  saveRDS(markers1, file = "./outdata/GSE149614/other/markers1.rds")

#### 2.2.4.1  abundance & malignant --------------------------------------------------------------
  ##single cell C1(1,2) C2(6,9,10,14)
  ##ssgsea  
  geneset <- split(markers1$gene, markers1$cluster)
  write.xlsx(geneset, file = "./outdata/GSE149614/geneset_can.xlsx")
  
  eggtest <- gsva(egg_us(Bulk_exp$TCGA),geneset,method='ssgsea')
  tumordata <- data.frame(predict_data$TCGA$Subtype,t(eggtest))
  colnames(tumordata) <- c("subtype", "HCC1","HCC2","HCC6","HCC9","HCC10","HCC14")
  plot_data <- reshape2::melt(tumordata, id.vars="subtype")
  plot_data$single_label <- ifelse(plot_data$variable %in% c("HCC6","HCC9","HCC10","HCC14"), "C1", "C2")
  #abundance
  p1 <- ggplot(plot_data,aes(x=interaction(variable,single_label),value,fill=subtype))+
    geom_violin(scale = "width",color=NA,position = position_dodge(1))+
    geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(1),
                 aes(color = subtype),fill="white")+
    scale_fill_manual(values=c("C1"="#CD1818","C2"="#0F2C67"))+
    scale_colour_manual(values=c("C1"="#CD1818","C2"="#0F2C67"))+
    theme_bw()+#theme(panel.grid=element_blank())+
    font("xy", size = 15)+
    font("xy.text", size = 14,color = "black")+
    font("legend.text", size = 15)+
    font("legend.title", size = 15)+
    theme(axis.title.x =element_blank())+
    stat_compare_means(aes(group = subtype),label = "p.signif",vjust=0.8,hide.ns=T,size=6
    )+
    scale_x_discrete(position = "top")+
    guides(x="axis_nested")+
    theme(strip.placement = "outside")+
    theme(strip.background.y = element_rect(fill = "grey85", colour = "grey85"),
          panel.border=element_rect(colour ="grey70"),
          strip.text.y = element_text(size = 14))+
    ylab("Abundance")
  
  ##correlation
  geneset <- freadlist("./data/cor_geneset.txt")
  library(GSVA)
  library(Hmisc)
  eggtest <- gsva(egg_us(Bulk_exp$TCGA),geneset,method='ssgsea')
  cor_data <- data.frame(t(eggtest),tumordata[,-1])
  cor_result  <-rcorr(as.matrix(cor_data),type = "spearman")
  cor_list <- list()
  cor_list$r <- cor_result$r[1:5,6:11]
  cor_list$p <- cor_result$P[1:5,6:11]
  
  corda <- data.frame((cor_list$p < 0.05) * cor_list$r,check.names=FALSE)
  corda$y <- rownames(corda)
  da <- reshape2::melt(data = corda) %>% na.omit()
  da[da$value==0,3] <- NA
  da$class <- ifelse(da$variable %in% c("HCC6","HCC9","HCC10","HCC14"),"C1","C2")
  da$variable <- factor(da$variable,levels = c("HCC1","HCC2","HCC6","HCC9","HCC10","HCC14"))
  da$class <- factor(da$class,levels = c("C1","C2"))
  da <- na.omit(da)
  
  p2 <- ggplot(da, aes(x=interaction(variable,class),y = y)) +
    geom_tile(aes(fill = value),color = 'white')+
    scale_fill_gradient2(low ="#0F2C67",high = "#CD1818",mid = "white",midpoint = 0,na.value = "white")+
    geom_text(aes(interaction(variable,class), y, label = round(value,2)), color = "black", size = 5)+
    guides(x="axis_nested")+   
    theme_bw()+
    theme(panel.grid = element_blank())+ 
    font("legend.text", size = 15)+
    font("legend.title", size = 15)+
    theme(
      panel.border=element_rect(colour ="white"),
      strip.text.y = element_text(size = 14),
      axis.title.x =element_blank(),
      axis.title.y =element_blank(),
      axis.text = element_text(color = 'black',size = 14),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x =  element_blank())
  
  #拼图
  library(aplot)
  pdf(file = "./plot/cancer/cor_bulk_score.pdf",width =6.5 ,height = 5)
  p1 %>% insert_bottom(p2)
  dev.off()
#### 2.2.4.2  survival --------------------------------------------------------------  
  library(survival)
  library(survminer)
  geneset <- split(markers1$gene, markers1$cluster)
  eggtest <- gsva(egg_us(Bulk_exp$TCGA),geneset,method='ssgsea')
  Bulk_data <- readRDS("./data/Bulk_data.rds")
  name = c("HCC1","HCC2","HCC6","HCC9","HCC10","HCC14")
  sur_plot <- list()
  for (i in 1:6) {
    sur_plot[i] <- feature_surplot(data = eggtest[i,], pdata = Bulk_data$TCGA$pdata,
                  time = "OS.time", event ="OS",count = F,
                  cutoff_method = "best", title = name[i])
  }
  
  pdf(file = "./plot/cancer/survival.pdf",width =11 ,height = 7)
  wrap_plots(sur_plot)
  dev.off()
  
  sur_HR <- list()
  for (i in 1:6) {
    sur_HR[[i]] <- feature_surplot(data = eggtest[i,], pdata = Bulk_data$TCGA$pdata,
                                   time = "OS.time", event ="OS",count = F,
                                   cutoff_method = "best", title = name[i],P_HR=T)
  }
  sur_HR <- do.call(rbind,sur_HR)
  sur_HR <- as.data.frame(sur_HR)
  sur_HR$variable <-name
  sur_HR$class <- ifelse(sur_HR$variable %in% c("HCC6","HCC9","HCC10","HCC14"),"C1","C2")
  sur_HR$y <- "HR"
  sur_HR$label <- paste0(round(sur_HR$HR,2), ifelse(sur_HR$p.val < 0.05 ,"*",""))
  
  p3 <- ggplot(sur_HR, aes(x=interaction(variable,class),y = y)) +
    geom_tile(aes(fill = HR),color = 'white')+
    scale_fill_gradient2(low ="#0F2C67",high = "#CD1818",mid = "white",midpoint = 1,na.value = "white")+
    geom_text(aes(interaction(variable,class), y, label = label), color = "black", size = 5)+
    guides(x="axis_nested")+   
    theme_bw()+
    theme(panel.grid = element_blank())+ 
    font("legend.text", size = 15)+
    font("legend.title", size = 15)+
    theme(
      panel.border=element_rect(colour ="white"),
      strip.text.y = element_text(size = 14),
      axis.title.x =element_blank(),
      axis.title.y =element_blank(),
      axis.text = element_text(color = 'black',size = 14),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x =  element_blank())
  
  pdf(file = "./plot/cancer/cor_bulk_score1.pdf",width =6.5 ,height = 5)
  p1 %>% insert_bottom(p2) %>% insert_bottom(p3,height=.25) 
  dev.off()
  
#### 2.2.4.3  cancer gene --------------------------------------------------------------  
  #####查看癌基因
  DTG<-fread("./data/NCG_cancerdrivers_systemslevelproperties.tsv")[,1:2]
  TSG<-fread("./data/TSgene.csv")[,1:2]
  DTG_marker<-markers1[match(intersect(markers1$gene,DTG$symbol),markers1$gene),]
  TSG_marker<-markers1[match(intersect(markers1$gene,TSG$GeneSymbol),markers1$gene),]
  table(DTG_marker$cluster)
  table(TSG_marker$cluster)
  ##画柱状图
  df<-rbind(data.frame(table(DTG_marker$cluster),type="DTG"),
            data.frame(table(TSG_marker$cluster),type="TSG"))
  df$Var1 <- paste0("HCC",df$Var1)
  pdf("./plot/cancer/DTG_TSG.pdf",width=3.5,height=3.5)
  ggbarplot(df, x="Var1", y="Freq",
            fill = "type", color = "type",
            label = TRUE, lab.col = "white", lab.pos = "in",
            palette=c("nejm"),
            width=0.8,
            ggtheme=theme_bw())+
    font("xy.text",size=14)+
    font("legend.text",size=14)+
    font("legend.title",size=14)+
    font("xy",size=14)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    labs(x="",y="Number")
  dev.off()
  

## 2.3 GSE151530 Tumor---------------------------------------------------------------
  
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
  source("./code/validation.R")
  Tumor <- readRDS("./outdata/GSE151530/Tumor.rds")
  Tumor <- NormalizeData(Tumor)
  data<-Tumor_pair(Tumor,result,result3)
  data<-na.omit(data)
  subtype <- Tumor_class(Tumor,data,group.by = "Sample")
  
  ##肿瘤细胞计算得分分类并画图
  Tumor <- tumor_score(Tumor,subtype = subtype, group.by="Sample")
  name=c("ISG.RS", "Proliferation","CTA","Stem","Hypoxia","Glycloysis",
         "SPP1","ITGA5","ITGB1")
  Tumormeta <- Tumor@meta.data[,c("subtype","Type","Proliferation","CTA","Stem","Hypoxia","Glycloysis")]
  Tumormeta <- reshape2::melt(Tumormeta,id.vars=c("subtype","Type"))
  #Tumormeta$seurat_clusters <- factor(Tumormeta$seurat_clusters, levels = c("HCC1","HCC2","HCC6","HCC9","HCC10","HCC14"))
  rm(Tumor)
  
  ##画图
  library(dittoSeq)
  library(patchwork)
  name <- as.character(unique(Tumormeta$variable))
  
  p1 <- plot.violin(Tumormeta,x=variable,y=value,
              fill = subtype,
              palette.fill=c("#CD1818","#0F2C67"),
              palette.color = c("#CD1818","#0F2C67"),
              color = "white",
              fill.b = "white",
              color.b = subtype,
              group=subtype,
              alpha.v = 0.7,
              notch.b=F,
              width.b=0.5,
              add = c("boxplot","violinplot","shadow"),#,"shadow"
              shadow.col = c("grey40","white"),
              lwd.b = 0.8,
              x_label="",
              y_label = "Score",
              theme = theme_bw(base_size = 16))+
    theme(
          plot.margin = margin(t = 1,unit = 'cm'))+
    theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  pdf("./plot/GSE151530_Tumor.pdf", width = 5, height = 4)
  p1
  dev.off()
  
  
# 3. Lymphocyte ---------------------------------------------------------------
## 3.1 Lymphocyte score ---------------------------------------------------------------
### 3.1.1 GSE149614 score ---------------------------------------------------------------
  # 处理淋巴 
  subtype <- readRDS("./pair/GSE149614_subtype.rds") 
  Lymphocyte4 <- readRDS("./outdata/GSE149614/lymphocyte/Lymphocyte4.rds")
  #Ly <- subset(Lymphocyte4, main %in% c("CD4","CD8","Treg"))
  Ly <- Lymphocyte4
  test <- ly_score(Ly, subtype, group.by = "Sample")
#### 3.1.1.1 Exhausted ---------------------------------------------------------------
  ##Exhausted
  library(ggh4x)
  test1 <- subset(test, main %in% c("CD4","CD8","Treg"))
  plot_data <- test1@meta.data[,c("Exhausted", "subtype", "main")]
  plot_data$type <- "Exhausted"
  plot_data$dataset <- "GSE149614"
  plot_data1 <- plot_data[plot_data$main %in% "CD8", ]
  plot_data$main <- "T cell"
  plot_data <- rbind(plot_data, plot_data1)
  plot_data$main <- factor(plot_data$main, levels = c("CD8","T cell"), labels = c("CD8+T","T cell" ))
  
  ##各亚型
  p1 <- plot.violin(plot_data,x=main,y=Exhausted,
              fill = subtype,
              palette.fill=c("#CD1818","#0F2C67"),
              palette.color = c("#CD1818","#0F2C67"),
              color = "white",
              fill.b = "white",
              color.b = subtype,
              group=subtype,
              alpha.v = 0.7,
              notch.b=F,
              width.b=0.4,
              add = c("boxplot","violinplot","shadow"),#,"shadow"
              lwd.b = 0.8,
              shadow.col = c("grey40","white"),
              x_label="",
              y_label = "",
              theme = theme_bw(base_size = 16))+
    theme(plot.margin = margin(t = 1,unit = 'cm'))+
    coord_cartesian(clip = 'off')+
    facet_wrap2("dataset",scales ="free_y")+
    labs(y="Exhausted",color = "Subtype",fill= "Subtype")
  
#### 3.1.1.2 costimulatory ---------------------------------------------------------------
  library(ggh4x)
  test1 <- subset(test, main %in% c("CD4","CD8","Treg"))
  plot_data <- test1@meta.data[,c("costimulatory", "subtype", "main")]
  plot_data$type <- "Costimulatory"
  plot_data$dataset <- "GSE149614"
  
  ##各亚型
  p2 <- plot.violin(plot_data,x=main,y=costimulatory,
                    fill = subtype,
                    palette.fill=c("#CD1818","#0F2C67"),
                    palette.color = c("#CD1818","#0F2C67"),
                    color = "white",
                    fill.b = "white",
                    color.b = subtype,
                    group=subtype,
                    alpha.v = 0.7,
                    notch.b=F,
                    width.b=0.4,
                    add = c("boxplot","violinplot","shadow"),#,"shadow"
                    lwd.b = 0.8,
                    shadow.col = c("grey40","white"),
                    x_label="",
                    y_label = "",
                    theme = theme_bw(base_size = 16))+
    theme(plot.margin = margin(t = 1,unit = 'cm'))+
    coord_cartesian(clip = 'off')+
    facet_wrap2("dataset",scales ="free_y")+
    labs(y="Costimulatory",color = "Subtype",fill= "Subtype")
  
### 3.1.2 GSE151530 Lymphocyte ---------------------------------------------------------------
  # 处理淋巴 
  subtype <- readRDS("./pair/GSE151530_subtype.rds") 
  Ly3 <- readRDS("./outdata/GSE151530/Ly3.rds")
  #Ly <- subset(Lymphocyte4, main %in% c("CD4","CD8","Treg"))
  Ly <- Ly3
  test <- ly_score(Ly, subtype, group.by = "Sample")
  
#### 3.1.2.1 Exhausted ---------------------------------------------------------------
  ##Exhausted
  library(ggh4x)
  test1 <- subset(test, main %in% c("CD4","CD8","Treg"))
  plot_data <- test1@meta.data[,c("Exhausted", "subtype", "main")]
  plot_data$type <- "Exhausted"
  plot_data$dataset <- "GSE151530"
  plot_data1 <- plot_data[plot_data$main %in% "CD8", ]
  plot_data$main <- "T cell"
  plot_data <- rbind(plot_data, plot_data1)
  plot_data$main <- factor(plot_data$main, levels = c("CD8","T cell" ), labels = c("CD8+T","T cell" ) )
  
  ##各亚型
  p3 <- plot.violin(plot_data,x=main,y=Exhausted,
                    fill = subtype,
                    palette.fill=c("#CD1818","#0F2C67"),
                    palette.color = c("#CD1818","#0F2C67"),
                    color = "white",
                    fill.b = "white",
                    color.b = subtype,
                    group=subtype,
                    alpha.v = 0.7,
                    notch.b=F,
                    width.b=0.4,
                    add = c("boxplot","violinplot","shadow"),#,"shadow"
                    lwd.b = 0.8,
                    shadow.col = c("grey40","white"),
                    x_label="",
                    y_label = "",
                    theme = theme_bw(base_size = 16))+
    theme(plot.margin = margin(t = 1,unit = 'cm'))+
    coord_cartesian(clip = 'off')+
    facet_wrap2("dataset",scales ="free_y")+
    labs(y="Exhausted",color = "Subtype",fill= "Subtype")
  
#### 3.1.2.2 costimulatory ---------------------------------------------------------------
  library(ggh4x)
  test1 <- subset(test, main %in% c("CD4","CD8","Treg"))
  plot_data <- test1@meta.data[,c("costimulatory", "subtype", "main")]
  plot_data$type <- "Costimulatory"
  plot_data$dataset <- "GSE151530"
  
  
  ##各亚型
  p4 <- plot.violin(plot_data,x=main,y=costimulatory,
                    fill = subtype,
                    palette.fill=c("#CD1818","#0F2C67"),
                    palette.color = c("#CD1818","#0F2C67"),
                    color = "white",
                    fill.b = "white",
                    color.b = subtype,
                    group=subtype,
                    alpha.v = 0.7,
                    notch.b=F,
                    width.b=0.4,
                    add = c("boxplot","violinplot","shadow"),#,"shadow"
                    lwd.b = 0.8,
                    shadow.col = c("grey40","white"),
                    x_label="",
                    y_label = "",
                    theme = theme_bw(base_size = 16))+
    theme(plot.margin = margin(t = 1,unit = 'cm'))+
    coord_cartesian(clip = 'off')+
    facet_wrap2("dataset",scales ="free_y")+
    labs(y="Costimulatory",color = "Subtype",fill= "Subtype")
  
### 3.1.3 save plot ---------------------------------------------------------------
  pdf("./plot/lymphocyte/score.pdf", width = 7, height = 5.5)
  (p1 + p2 + plot_layout(guides = "collect",widths = c(2,3)))/
  (p3 + p4 + plot_layout(guides = "collect",widths = c(2,3)))
  dev.off()
  
## 3.2  Marker genes ( abundance, malignant, survival) ---------------------------------------------------------------
  library(GSVA)
  library(ggpubr)
  library(ggh4x)
  Bulk_exp <- readRDS("./data/Bulk_exp.rds")
  predict_data <- readRDS("./data/predict_data.rds")
  ##marker gene
  # 处理淋巴 
  subtype <- readRDS("./pair/GSE149614_subtype.rds") 
  Lymphocyte4 <- readRDS("./outdata/GSE149614/lymphocyte/Lymphocyte4.rds")
  #Ly <- subset(Lymphocyte4, main %in% c("CD4","CD8","Treg"))
  Ly <- Lymphocyte4
  test <- ly_score(Ly, subtype, group.by = "Sample")
  #test1 <- subset(test, main %in% c("CD4","CD8","Treg"))
  marker<-FindAllMarkers(test, only.pos = T,min.pct=0.4,logfc.threshold = 0.5)
  markers1 <- marker[!grepl("^RP[ SL ]", marker$gene, ignore.case = F),]
  markers1 <- markers1[!grepl("^MT-", markers1$gene, ignore.case = F),]
  #saveRDS(markers1, file = "./outdata/GSE149614/other/markers1.rds")
  markers2 <- markers1[markers1$cluster %in% c("0","1","2","4A","4B","6A","6B"), ]
  markers2$cluster <- as.character(markers2$cluster)
  
### 3.2.1  abundance & malignant --------------------------------------------------------------
  ##ssgsea  
  geneset <- split(markers2$gene, markers2$cluster)
  library(openxlsx)
  write.xlsx(geneset, file = "./outdata/GSE149614/geneset_ly.xlsx")
  eggtest <- gsva(egg_us(Bulk_exp$TCGA),geneset,method='ssgsea')
  data <- data.frame(predict_data$TCGA$Subtype,t(eggtest))
  colnames(data) <- c("subtype", paste0("L",rownames(eggtest)))
  plot_data <- reshape2::melt(data, id.vars="subtype")
  plot_data$single_label <- ifelse(plot_data$variable %in% c("L0"), "C2", "C1")
  #abundance
  p1 <- ggplot(plot_data,aes(x=interaction(variable,single_label),value,fill=subtype))+
    geom_violin(scale = "width",color=NA,position = position_dodge(1))+
    geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(1),
                 aes(color = subtype),fill="white")+
    scale_fill_manual(values=c("C1"="#CD1818","C2"="#0F2C67"))+
    scale_colour_manual(values=c("C1"="#CD1818","C2"="#0F2C67"))+
    theme_bw()+#theme(panel.grid=element_blank())+
    font("xy", size = 15)+
    font("xy.text", size = 14,color = "black")+
    font("legend.text", size = 15)+
    font("legend.title", size = 15)+
    theme(axis.title.x =element_blank())+
    stat_compare_means(aes(group = subtype),label = "p.signif",vjust=0.8,hide.ns=T,size=6
    )+
    scale_x_discrete(position = "top")+
    guides(x="axis_nested")+
    theme(strip.placement = "outside")+
    theme(strip.background.y = element_rect(fill = "grey85", colour = "grey85"),
          panel.border=element_rect(colour ="grey70"),
          strip.text.y = element_text(size = 14))+
    ylab("Abundance")
  p1
  ##correlation
  geneset <- freadlist("./data/cor_geneset.txt")
  library(GSVA)
  library(Hmisc)
  eggtest <- gsva(egg_us(Bulk_exp$TCGA),geneset,method='ssgsea')
  cor_data <- data.frame(t(eggtest),data[,-1])
  cor_result  <-rcorr(as.matrix(cor_data),type = "spearman")
  cor_list <- list()
  cor_list$r <- cor_result$r[1:5,6:12]
  cor_list$p <- cor_result$P[1:5,6:12]
  
  corda <- data.frame((cor_list$p < 0.05) * cor_list$r,check.names=FALSE)
  corda$y <- rownames(corda)
  da <- reshape2::melt(data = corda) %>% na.omit()
  da[da$value==0,3] <- NA
  da$class <- ifelse(da$variable %in% c("L0"),"C2","C1")
  #da$variable <- factor(da$variable,levels = c("HCC1","HCC2","HCC6","HCC9","HCC10","HCC14"))
  da$class <- factor(da$class,levels = c("C1","C2"))
  da <- na.omit(da)
  da$y <- factor(da$y, levels = c("Glycloysis", "Hypoxia", "Proliferation","Stem","CTA" )[5:1])
  
  p2 <- ggplot(da, aes(x=interaction(variable,class),y = y)) +
    geom_tile(aes(fill = value),color = 'white')+
    scale_fill_gradient2(low ="#0F2C67",high = "#CD1818",mid = "white",midpoint = 0,na.value = "white")+
    geom_text(aes(interaction(variable,class), y, label = round(value,2)), color = "black", size = 5)+
    guides(x="axis_nested")+   
    theme_bw()+
    theme(panel.grid = element_blank())+ 
    font("legend.text", size = 15)+
    font("legend.title", size = 15)+
    theme(
      panel.border=element_rect(colour ="white"),
      strip.text.y = element_text(size = 14),
      axis.title.x =element_blank(),
      axis.title.y =element_blank(),
      axis.text = element_text(color = 'black',size = 14),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x =  element_blank())
  
  #拼图
  library(aplot)
  pdf(file = "./plot/lymphocyte/cor_bulk_score.pdf",width =6.5 ,height = 4.5)
  p1 %>% insert_bottom(p2)
  dev.off()
### 3.2.2  survival --------------------------------------------------------------  
  library(survival)
  library(survminer)
  geneset <- split(markers2$gene, markers2$cluster)
  eggtest <- gsva(egg_us(Bulk_exp$TCGA),geneset,method='ssgsea')
  Bulk_data <- readRDS("./data/Bulk_data.rds")
  name = c(paste0("L",rownames(eggtest)))
  rownames(eggtest)<-name
  sur_plot <- list()
  for (i in 1:7) {
    sur_plot[i] <- feature_surplot(data = eggtest[i,], pdata = Bulk_data$TCGA$pdata,
                                   time = "OS.time", event ="OS",count = F,
                                   cutoff_method = "best", title = name[i])
  }
  
  pdf(file = "./plot/lymphocyte/survival.pdf",width =11 ,height = 10)
  wrap_plots(sur_plot)
  dev.off()
  
  sur_HR <- list()
  for (i in 1:7) {
    sur_HR[[i]] <- feature_surplot(data = eggtest[i,], pdata = Bulk_data$TCGA$pdata,
                                   time = "OS.time", event ="OS",count = F,
                                   cutoff_method = "best", title = name[i],P_HR=T)
  }
  sur_HR <- do.call(rbind,sur_HR)
  sur_HR <- as.data.frame(sur_HR)
  sur_HR$variable <-name
  sur_HR$class <- ifelse(sur_HR$variable %in% c("L0"),"C2","C1")
  sur_HR$y <- "HR"
  sur_HR$label <- paste0(round(sur_HR$HR,2), ifelse(sur_HR$p.val < 0.05 ,"*",""))
  
  p3 <- ggplot(sur_HR, aes(x=interaction(variable,class),y = y)) +
    geom_tile(aes(fill = HR),color = 'white')+
    scale_fill_gradient2(low ="#0F2C67",high = "#CD1818",mid = "white",midpoint = 1,na.value = "white")+
    geom_text(aes(interaction(variable,class), y, label = label), color = "black", size = 5)+
    guides(x="axis_nested")+   
    theme_bw()+
    theme(panel.grid = element_blank())+ 
    font("legend.text", size = 15)+
    font("legend.title", size = 15)+
    theme(
      panel.border=element_rect(colour ="white"),
      strip.text.y = element_text(size = 14),
      axis.title.x =element_blank(),
      axis.title.y =element_blank(),
      axis.text = element_text(color = 'black',size = 14),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x =  element_blank())
  
  pdf(file = "./plot/lymphocyte/cor_bulk_score1.pdf",width =6.5 ,height = 5)
  p1 %>% insert_bottom(p2) %>% insert_bottom(p3,height=.25) 
  dev.off()
# 4. Myeloid ---------------------------------------------------------------
## 4.1 Myeloid score ---------------------------------------------------------------
### 4.1.1 GSE149614 score ---------------------------------------------------------------
  # 处理髓系
  subtype <- readRDS("./pair/GSE149614_subtype.rds") 
  My <- readRDS("./outdata/GSE149614/myeloid/Myeloid4.rds")
  test <- my_score(My, subtype, group.by = "Sample")
#### 4.1.1.1 GSE149614 Myeloid ---------------------------------------------------------------
  library(dplyr)
  
  test$ano <- paste0("M",test$seurat_clusters)
  test$ano <- factor(test$ano, levels = c(paste0("M",0:11)))
  test1 <- subset(test, main %in% c("Macrophage"))
  test1$dataset <- "GSE149614"
  ##EMR
  a<-dittoPlot(test1, "EMR",group.by = "ano",color.panel = dittoColors()[c(0,2,4,8,1,6)+1] ,
            plots = c("boxplot"),
            boxplot.width = 0.5)
  a1 <- a$data %>% group_by(grouping) %>% summarise(median=median(var.data), mean=mean(var.data))
  p1 <- dittoPlot(test1, "EMR",group.by = "ano",color.panel = dittoColors()[c(0,2,4,8,1,6)+1] ,
            plots = c("boxplot"),x.reorder=order(a1$median, a1$mean,decreasing = T),
            boxplot.width = 0.7,theme=theme_bw(),xlab=NULL,main ="EMR GSE149614")
  
  ##Proangiogenic
  a<-dittoPlot(test1, "Proangiogenic",group.by = "ano",color.panel = dittoColors()[c(0,2,4,8,1,6)+1] ,
               plots = c("boxplot"),
               boxplot.width = 0.5)
  a1 <- a$data %>% group_by(grouping) %>% summarise(median=median(var.data), mean=mean(var.data))
  p2 <- dittoPlot(test1, "Proangiogenic",group.by = "ano",color.panel = dittoColors()[c(0,2,4,8,1,6)+1] ,
            plots = c("boxplot"),x.reorder=order(a1$median, a1$mean,decreasing = T),
            boxplot.width = 0.7,theme=theme_bw(),xlab=NULL,main ="Proangiogenic GSE149614")
  
  
  ##Lipid
  a<-dittoPlot(test1, "Lipid",group.by = "ano",color.panel = dittoColors()[c(0,2,4,8,1,6)+1] ,
               plots = c("boxplot"),
               boxplot.width = 0.5)
  a1 <- a$data %>% group_by(grouping) %>% summarise(median=median(var.data), mean=mean(var.data))
  p3<-dittoPlot(test1, "Lipid",group.by = "ano",color.panel = dittoColors()[c(0,2,4,8,1,6)+1] ,
            plots = c("boxplot"),x.reorder=order(a1$median, a1$mean,decreasing = T),
            boxplot.width = 0.7,theme=theme_bw(),xlab=NULL,main ="Lipid GSE149614")
  
  #MHC
  a<-dittoPlot(test, "MHC2",group.by = "ano",#color.panel = dittoColors()[c(0,2,4,8,1,6)+1] ,
               plots = c("boxplot"),
               boxplot.width = 0.5)
  a1 <- a$data %>% group_by(grouping) %>% summarise(median=median(var.data), mean=mean(var.data))
  p4<-dittoPlot(test, "MHC2",group.by = "ano",#color.panel = dittoColors()[c(0,2,4,8,1,6)+1] ,
            plots = c("boxplot"),x.reorder=order(a1$median, a1$mean,decreasing = T),
            boxplot.width = 0.7,theme=theme_bw(),xlab=NULL,main ="MHC2 GSE149614")
  
  
  
### 4.1.2 GSE151530 score ---------------------------------------------------------------
  # 处理髓系
  subtype <- readRDS("./pair/GSE151530_subtype.rds") 
  My <- readRDS("./outdata/GSE151530/Myeloid2.rds")
  test <- my_score(My, subtype, group.by = "Sample")
#### 4.1.2.1 GSE151530 Myeloid ---------------------------------------------------------------
  library(dplyr)
  
  test$ano <- paste0("M",test$seurat_clusters)
  test$ano <- factor(test$ano, levels = c(paste0("M",0:14)))
  test1 <- subset(test, main %in% c("Macrophage"))
  
  ##EMR
  a<-dittoPlot(test1, "EMR",group.by = "ano",#color.panel = dittoColors()[c(0,2,4,8,1,6)+1] ,
               plots = c("boxplot"),
               boxplot.width = 0.5)
  a1 <- a$data %>% group_by(grouping) %>% summarise(median=median(var.data), mean=mean(var.data))
  p5 <- dittoPlot(test1, "EMR",group.by = "ano",#color.panel = dittoColors()[c(0,2,4,8,1,6)+1] ,
            plots = c("boxplot"),x.reorder=order(a1$median, a1$mean,decreasing = T),
            boxplot.width = 0.7,theme=theme_bw(),xlab=NULL,main ="EMR GSE151530")
  
  ##Proangiogenic
  a<-dittoPlot(test1, "Proangiogenic",group.by = "ano",#color.panel = dittoColors()[c(0,2,4,8,1,6)+1] ,
               plots = c("boxplot"),
               boxplot.width = 0.5)
  a1 <- a$data %>% group_by(grouping) %>% summarise(median=median(var.data), mean=mean(var.data))
  p6 <- dittoPlot(test1, "Proangiogenic",group.by = "ano",#color.panel = dittoColors()[c(0,2,4,8,1,6)+1] ,
            plots = c("boxplot"),x.reorder=order(a1$median, a1$mean,decreasing = T),
            boxplot.width = 0.7,theme=theme_bw(),xlab=NULL,main ="Proangiogenic GSE151530")
  
  
  ##Lipid
  a<-dittoPlot(test1, "Lipid",group.by = "ano",#color.panel = dittoColors()[c(0,2,4,8,1,6)+1] ,
               plots = c("boxplot"),
               boxplot.width = 0.5)
  a1 <- a$data %>% group_by(grouping) %>% summarise(median=median(var.data), mean=mean(var.data))
  p7 <- dittoPlot(test1, "Lipid",group.by = "ano",#color.panel = dittoColors()[c(0,2,4,8,1,6)+1] ,
            plots = c("boxplot"),x.reorder=order(a1$median, a1$mean,decreasing = T),
            boxplot.width = 0.7,theme=theme_bw(),xlab=NULL,main ="Lipid GSE151530")
  
  #MHC2
  a<-dittoPlot(test, "MHC2",group.by = "ano",#color.panel = dittoColors()[c(0,2,4,8,1,6)+1] ,
               plots = c("boxplot"),
               boxplot.width = 0.5)
  a1 <- a$data %>% group_by(grouping) %>% summarise(median=median(var.data), mean=mean(var.data))
  p8 <- dittoPlot(test, "MHC2",group.by = "ano",#color.panel = dittoColors()[c(0,2,4,8,1,6)+1] ,
            plots = c("boxplot"),x.reorder=order(a1$median, a1$mean,decreasing = T),
            boxplot.width = 0.7,theme=theme_bw(),xlab=NULL,main ="MHC2 GSE151530")
  
### 4.1.3 patch plot ---------------------------------------------------------------
  library(patchwork)
  pdf("./plot/myeloid/MA_score.pdf", width = 8, height = 4.5)
  (p1 + p3 + p2 + plot_layout(guides = "collect"))/
  (p5 + p7 + p6 + plot_layout(guides = "collect"))
  dev.off()
  
  pdf("./plot/myeloid/MHC2.pdf", width = 5, height = 5)
  p4/p8
  dev.off()
  
## 4.2  Marker genes ( abundance, malignant, survival) ---------------------------------------------------------------
  library(GSVA)
  library(ggpubr)
  library(ggh4x)
  Bulk_exp <- readRDS("./data/Bulk_exp.rds")
  predict_data <- readRDS("./data/predict_data.rds")
  ##marker gene
  # 处理髓系
  subtype <- readRDS("./pair/GSE149614_subtype.rds") 
  My <- readRDS("./outdata/GSE149614/myeloid/Myeloid4.rds")
  test <- my_score(My, subtype, group.by = "Sample")
  test1 <- subset(test, main %in% c("Macrophage"))
  marker<-FindAllMarkers(test, only.pos = T,min.pct = 0.2,logfc.threshold = 0.5)
  markers1 <- marker[!grepl("^RP[ SL ]", marker$gene, ignore.case = F),]
  markers1 <- markers1[!grepl("^MT-", markers1$gene, ignore.case = F),]
  #saveRDS(markers1, file = "./outdata/GSE149614/other/markers1.rds")
  markers2 <- markers1[markers1$cluster %in% c("0","2","4","8","1","6"), ]
  markers2$cluster <- as.character(markers2$cluster)
### 4.2.1  abundance & malignant --------------------------------------------------------------
  ##ssgsea  
  geneset <- split(markers2$gene, markers2$cluster)
  write.xlsx(geneset, file = "./outdata/GSE149614/geneset_my.xlsx")
  
  eggtest <- gsva(egg_us(Bulk_exp$TCGA),geneset,method='ssgsea')
  data <- data.frame(predict_data$TCGA$Subtype,t(eggtest))
  colnames(data) <- c("subtype", paste0("M",rownames(eggtest)))
  plot_data <- reshape2::melt(data, id.vars="subtype")
  plot_data$single_label <- ifelse(plot_data$variable %in% c("M1","M6"), "C1", "C2")
  #abundance
  p1 <- ggplot(plot_data,aes(x=interaction(variable,single_label),value,fill=subtype))+
    geom_violin(scale = "width",color=NA,position = position_dodge(1))+
    geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(1),
                 aes(color = subtype),fill="white")+
    scale_fill_manual(values=c("C1"="#CD1818","C2"="#0F2C67"))+
    scale_colour_manual(values=c("C1"="#CD1818","C2"="#0F2C67"))+
    theme_bw()+#theme(panel.grid=element_blank())+
    font("xy", size = 15)+
    font("xy.text", size = 14,color = "black")+
    font("legend.text", size = 15)+
    font("legend.title", size = 15)+
    theme(axis.title.x =element_blank())+
    stat_compare_means(aes(group = subtype),label = "p.signif",vjust=0.8,hide.ns=T,size=6
    )+
    scale_x_discrete(position = "top")+
    guides(x="axis_nested")+
    theme(strip.placement = "outside")+
    theme(strip.background.y = element_rect(fill = "grey85", colour = "grey85"),
          panel.border=element_rect(colour ="grey70"),
          strip.text.y = element_text(size = 14))+
    ylab("Abundance")
  p1
  ##correlation
  geneset <- freadlist("./data/cor_geneset.txt")
  library(GSVA)
  library(Hmisc)
  eggtest <- gsva(egg_us(Bulk_exp$TCGA),geneset,method='ssgsea')
  cor_data <- data.frame(t(eggtest),data[,-1])
  cor_result  <-rcorr(as.matrix(cor_data),type = "spearman")
  cor_list <- list()
  cor_list$r <- cor_result$r[1:5,6:11]
  cor_list$p <- cor_result$P[1:5,6:11]
  
  corda <- data.frame((cor_list$p < 0.05) * cor_list$r,check.names=FALSE)
  corda$y <- rownames(corda)
  da <- reshape2::melt(data = corda) %>% na.omit()
  da[da$value==0,3] <- NA
  da$class <- ifelse(da$variable %in% c("M1","M6"),"C1","C2")
  da$class <- factor(da$class,levels = c("C1","C2"))
  da <- na.omit(da)
  da$y <- factor(da$y, levels = c("Glycloysis", "Hypoxia", "Proliferation","Stem","CTA" )[5:1])
  
  p2 <- ggplot(da, aes(x=interaction(variable,class),y = y)) +
    geom_tile(aes(fill = value),color = 'white')+
    scale_fill_gradient2(low ="#0F2C67",high = "#CD1818",mid = "white",midpoint = 0,na.value = "white")+
    geom_text(aes(interaction(variable,class), y, label = round(value,2)), color = "black", size = 5)+
    guides(x="axis_nested")+   
    theme_bw()+
    theme(panel.grid = element_blank())+ 
    font("legend.text", size = 15)+
    font("legend.title", size = 15)+
    theme(
      panel.border=element_rect(colour ="white"),
      strip.text.y = element_text(size = 14),
      axis.title.x =element_blank(),
      axis.title.y =element_blank(),
      axis.text = element_text(color = 'black',size = 14),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x =  element_blank())
  
  #拼图
  library(aplot)
  pdf(file = "./plot/myeloid/cor_bulk_score.pdf",width =6.5 ,height = 4.5)
  p1 %>% insert_bottom(p2)
  dev.off()
  
### 4.2.2  survival --------------------------------------------------------------  
  library(survival)
  library(survminer)
  geneset <- split(markers2$gene, markers2$cluster)
  eggtest <- gsva(egg_us(Bulk_exp$TCGA),geneset,method='ssgsea')
  Bulk_data <- readRDS("./data/Bulk_data.rds")
  name = c(paste0("M",rownames(eggtest)))
  rownames(eggtest)<-name
  sur_plot <- list()
  for (i in 1:6) {
    sur_plot[i] <- feature_surplot(data = eggtest[i,], pdata = Bulk_data$TCGA$pdata,
                                   time = "OS.time", event ="OS",count = F,
                                   cutoff_method = "best", title = name[i])
  }
  
  pdf(file = "./plot/myeloid/survival.pdf",width =11 ,height = 7)
  wrap_plots(sur_plot)
  dev.off()  
  
  sur_HR <- list()
  for (i in 1:6) {
    sur_HR[[i]] <- feature_surplot(data = eggtest[i,], pdata = Bulk_data$TCGA$pdata,
                                   time = "OS.time", event ="OS",count = F,
                                   cutoff_method = "best", title = name[i],P_HR=T)
  }
  sur_HR <- do.call(rbind,sur_HR)
  sur_HR <- as.data.frame(sur_HR)
  sur_HR$variable <-name
  sur_HR$class <- ifelse(sur_HR$variable %in% c("M1","M6"),"C1","C2")
  sur_HR$y <- "HR"
  sur_HR$label <- paste0(round(sur_HR$HR,2), ifelse(sur_HR$p.val < 0.05 ,"*",""))
  
  p3 <- ggplot(sur_HR, aes(x=interaction(variable,class),y = y)) +
    geom_tile(aes(fill = HR),color = 'white')+
    scale_fill_gradient2(low ="#0F2C67",high = "#CD1818",mid = "white",midpoint = 1,na.value = "white")+
    geom_text(aes(interaction(variable,class), y, label = label), color = "black", size = 5)+
    guides(x="axis_nested")+   
    theme_bw()+
    theme(panel.grid = element_blank())+ 
    font("legend.text", size = 15)+
    font("legend.title", size = 15)+
    theme(
      panel.border=element_rect(colour ="white"),
      strip.text.y = element_text(size = 14),
      axis.title.x =element_blank(),
      axis.title.y =element_blank(),
      axis.text = element_text(color = 'black',size = 14),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x =  element_blank())
  
  pdf(file = "./plot/myeloid/cor_bulk_score.pdf",width =6.5 ,height = 5)
  p1 %>% insert_bottom(p2) %>% insert_bottom(p3,height=.25) 
  dev.off()
# 5. Correlations between different single-cell data-----------------------------------------------------------
## 5.1 lymphocyte --------------------------------------------------------------
  #lymphocyte 
  ##load data
  Ly_GSE149614 <- readRDS("E:/project/HCC_immune_scRNA/outdata/GSE149614/lymphocyte/Lymphocyte4.rds")
  Ly_GSE151530 <- readRDS("E:/project/HCC_immune_scRNA/outdata/GSE151530/Ly3.rds")
  Ly_list <- list()
  Ly_list$GSE149614 <- AverageExpression(Ly_GSE149614, group.by = "ano_cluster")$RNA
  Ly_list$GSE151530 <- AverageExpression(Ly_GSE151530, group.by = "seurat_clusters")$RNA
  colnames(Ly_list$GSE149614) <- paste0("GSE149614_L",colnames(Ly_list$GSE149614))
  colnames(Ly_list$GSE151530) <- paste0("GSE151530_L",colnames(Ly_list$GSE151530))
  
  #提取基因名取交集
  Ly_name <- lapply(Ly_list, function(x){
    rownames(x)
  })
  Ly_name1 <- Reduce(intersect,Ly_name)
  Ly_list <- lapply(Ly_list, function(x){
    x[match(Ly_name1,rownames(x)),]
  })
  
  
  ##GSE151530 与gse149614
  corr_matrix <- cor(cbind(Ly_list$GSE149614, Ly_list$GSE151530),method = "spearman")
  corr_matrix1 <- corr_matrix[1:13,14:22]
  corr_matrix1 <- t(corr_matrix[1:13,14:22])
  corr_matrix1<-as.data.frame(corr_matrix1)
  corr_matrix2 <- corr_matrix1[c("GSE151530_L0","GSE151530_L5","GSE151530_L7","GSE151530_L1","GSE151530_L2","GSE151530_L8"),
                               c("GSE149614_L0","GSE149614_L1","GSE149614_L4A","GSE149614_L4B")]
  library(pheatmap)
  p1<-pheatmap::pheatmap(corr_matrix2,display_numbers=T,number_color="black",fontsize_number=13,cluster_cols = F,
                     cluster_rows = F,number_format="%.3f", border_color="white",
                     color = c(colorRampPalette(colors = c("#F5D0D0","#CD1818"))(20)))
  pdf("./plot/lymphocyte/cor_cell.pdf", width = 4,height = 5)
  print(p1)
  dev.off()
  
## 5.2 myeloid cells --------------------------------------------------------------------
  ##load data
  My_GSE149614 <- readRDS("E:/project/HCC_immune_scRNA/outdata/GSE149614/myeloid/Myeloid4.rds")
  My_GSE151530 <- readRDS("E:/project/HCC_immune_scRNA/outdata/GSE151530/Myeloid2.rds")
  My_list <- list()
  My_list$GSE149614 <- AverageExpression(My_GSE149614, group.by = "seurat_clusters")$RNA
  My_list$GSE151530 <- AverageExpression(My_GSE151530, group.by = "seurat_clusters")$RNA
  colnames(My_list$GSE149614) <- paste0("GSE149614_M",colnames(My_list$GSE149614))
  colnames(My_list$GSE151530) <- paste0("GSE151530_M",colnames(My_list$GSE151530))
  #提取基因名取交集
  My_name <- lapply(My_list, function(x){
    rownames(x)
  })
  My_name1 <- Reduce(intersect,My_name)
  My_list <- lapply(My_list, function(x){
    x[match(My_name1,rownames(x)),]
  })
  
  
  
  corr_matrix <- cor(cbind(My_list$GSE149614, My_list$GSE151530),method = "spearman")
  corr_matrix1 <- t(corr_matrix[1:12,13:26])
  #corr_matrix1 <- t(corr_matrix[1:12,13:22])
  corr_matrix1<-as.data.frame(corr_matrix1)
  corr_matrix2 <- corr_matrix1[c("GSE151530_M0","GSE151530_M3","GSE151530_M6","GSE151530_M9","GSE151530_M7","GSE151530_M11","GSE151530_M12"),
                               c("GSE149614_M1","GSE149614_M6")]
  
  library(pheatmap)
  p2 <- pheatmap::pheatmap(corr_matrix2,display_numbers=T,number_color="black",fontsize_number=13,cluster_cols = F,
                     cluster_rows = F,number_format="%.3f", border_color="white",
                     color = c(colorRampPalette(colors = c("#F5D0D0","#CD1818"))(20)))
  
  pdf("./plot/myeloid/cor_cell.pdf", width = 3,height = 5)
  print(p2)
  dev.off()
  


# 6.infercnv---------------------------------------------------------------------
  library(infercnv)
  library(devtools)
  library(AnnoProbe)
  ###initial data
  Tumor <- readRDS("./outdata/GSE149614/Tumor.rds")  
  Myeloid4 <- readRDS("./outdata/GSE149614/myeloid/Myeloid4.rds")
  
  dfcount<-cbind(Tumor@assays[["RNA"]]@counts,Myeloid4@assays[["RNA"]]@counts)
  # 第二个文件样本信息矩阵
  groupinfo= data.frame(cellId = colnames(dfcount))
  groupinfo$cellType=c(paste0("HCC",Tumor$seurat_clusters),Myeloid4$main)
  rownames(groupinfo) <- groupinfo$cellId
  groupinfo<-groupinfo[,-1,drop=F]
  ##第三个文件
  library(AnnoProbe)
  geneInfor=annoGene(rownames(dfcount),"SYMBOL",'human')
  geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
  geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
  rownames(geneInfor) <- geneInfor$SYMBOL
  geneInfor<-geneInfor[,-1,drop=F]
  
  ## 这里可以去除性染色体
  # 也可以把染色体排序方式改变
  dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
  dfcount =dfcount [match( geneInfor[,1], rownames(dfcount) ),]
  ##
  # 输出
  expFile='./inf/expFile.txt'
  write.table(dfcount ,file = expFile,sep = '\t',quote = F)
  groupFiles='./inf/groupFiles.txt'
  head(groupinfo)
  write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
  head(geneInfor)
  geneFile='./inf/geneFile.txt'
  write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dfcount,
                                      annotations_file=groupinfo,
                                      gene_order_file= geneInfor,
                                      ref_group_names=c("Macrophage","Monocyte", "DC", "Mast")) # 如果有正常细胞的话，把正常细胞的分组填进去
  rm(dfcount)
  saveRDS(infercnv_obj,file = "./inf/infercnv_obj.rds")
  load("./inf/infercnv_obj.Rdata")
  gc()
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, 
                               out_dir="./inf/inf3",
                               cluster_by_groups=TRUE, 
                               #analysis_mode="subclusters", #默认是"samples"
                               denoise=T,
                               HMM=F,output_format="pdf")
  
  
  #######cnv比例##
  library(dplyr)
  library(ggplot2)
  a <- fread(".\\inf\\inf3\\infercnv.preliminary.observations.txt",data.table = F)
  rownames(a) <- a[,1]
  a <- a[,-1]
  expr2=a-1
  expr2=expr2 ^ 2
  CNV_score=as.data.frame(colMeans(expr2))
  colnames(CNV_score)="CNV_score"
  CNV_score$CB=rownames(CNV_score)
  hmetadata <- Tumor@meta.data
  hmetadata$CB=rownames(hmetadata)
  CNV_score=CNV_score%>%inner_join(hmetadata,by="CB")
  
  ##main
  source("./code/Function.R")
  library(dittoSeq)
  library(ggpubr)
  CNV_score$seurat_clusters <- paste0("HCC", Tumor$seurat_clusters)
  CNV_score$seurat_clusters <- factor(CNV_score$seurat_clusters, levels = c("HCC1","HCC2","HCC6","HCC9","HCC10","HCC14"))
  CNV_score$subtype <- ifelse(CNV_score$Sample %in% c("HCC01T","HCC02T","HCC03T",
                                  "HCC04T","HCC05T","HCC06T"),"C2","C1")
  #cluster
  p1 <- plot.violin(CNV_score,x=seurat_clusters,y=CNV_score,
              fill = seurat_clusters,
              palette.fill=dittoColors()[c(1,2,6,9,10,14)+1],
              palette.color = dittoColors()[c(1,2,6,9,10,14)+1],
              color = "white",
              fill.b = "white",
              color.b = seurat_clusters,
              #group=seurat_clusters,
              alpha.v = 0.7,
              notch.b=F,
              width.b=0.5,
              lwd.b = 0.8,
              add = c("boxplot","violinplot"),#,"shadow"
              x_label="",
              y_label = "",
              theme = theme_bw(base_size = 16))+
    theme(axis.text.x = element_blank(),   ## 删去所有刻度标签
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 1,unit = 'cm'))+
    theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5))+
    coord_cartesian(clip = 'off')
  
  p2 <- plot.violin(CNV_score,x=seurat_clusters,y=CNV_score,
              fill = subtype,
              palette.fill=c("#CD1818","#0F2C67"),
              palette.color = c("#CD1818","#0F2C67"),
              color = "white",
              fill.b = "white",
              color.b = subtype,
              group=subtype,
              alpha.v = 0.7,
              notch.b=F,
              width.b=0.5,
              add = c("boxplot","violinplot","shadow"),#,"shadow"
              lwd.b = 0.8,
              shadow.col = c("grey40","white"),
              x_label="",
              y_label = "",
              theme = theme_bw(base_size = 16))+
    theme(#axis.text.x = element_blank(),   ## 删去所有刻度标签
          #axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),   ## 删去所有刻度标签
          axis.ticks.y = element_blank(),
          plot.margin = margin(t = 1,unit = 'cm'))+
    theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5))+
    coord_cartesian(clip = 'off')+
    theme(legend.position = 'none')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p3 <- plot.violin(CNV_score,x=Main,y=CNV_score,
              fill = subtype,
              palette.fill=c("#CD1818","#0F2C67"),
              palette.color = c("#CD1818","#0F2C67"),
              color = "white",
              fill.b = "white",
              color.b = subtype,
              group=subtype,
              alpha.v = 0.7,
              notch.b=F,
              width.b=0.5,
              add = c("boxplot","violinplot"),#,"shadow"
              lwd.b = 0.8,
              x_label="",
              y_label = "",
              theme = theme_bw(base_size = 16))+
    theme(axis.text.x = element_blank(),   ## 删去所有刻度标签
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),   ## 删去所有刻度标签
          axis.ticks.y = element_blank(),
          plot.margin = margin(t = 1,unit = 'cm'))+
    theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5))+
    coord_cartesian(clip = 'off')
  library(patchwork)
  
  pdf(file = "./plot/cancer/cnv.pdf",width =10 ,height = 3)
  p1+p2+p3+plot_layout(widths = c(2,4,1))
  dev.off()
  
  
  
  
  