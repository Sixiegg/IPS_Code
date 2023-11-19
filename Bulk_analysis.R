library(data.table)
source("./code/Function.R")
install.packages("./IPS29GPS_0.1.0.tar.gz", repos = NULL, type = "source")
# 1. Predict immuno-prognostic subtypes----------------------------------------------------------------------
library(IPS29GPS)
Bulk_exp <- readRDS("./data/Bulk_exp.rds")
pdata <- readRDS("./data/pdata.rds")
###predict
TCGA <- predict_29pairs(Bulk_exp$TCGA,Pairs)
ICGC <- predict_29pairs(Bulk_exp$ICGC,Pairs)
GSE14520 <- predict_29pairs(Bulk_exp$GSE14520,Pairs)
GSE116174 <- predict_29pairs(Bulk_exp$GSE116174,Pairs)
In_house <- predict_29pairs(Bulk_exp$"In-house",Pairs)
##
predict_data <- list(TCGA = TCGA,ICGC = ICGC,GSE14520 = GSE14520,GSE116174 = GSE116174,"In-house"=In_house)
Bulk_data <- list(TCGA=list(predict=TCGA, pdata=pdata$TCGA),
                  ICGC=list(predict=ICGC, pdata=pdata$ICGC),
                  GSE14520=list(predict=GSE14520, pdata=pdata$GSE14520),
                  GSE116174=list(predict=GSE116174, pdata=pdata$GSE116174),
                  "In-house"=list(predict=In_house, pdata=pdata$"In-house"))
##save
saveRDS(predict_data, file = "./outdata/predict_data.rds")
saveRDS(Bulk_data, file = "./outdata/Bulk_data.rds")


# 2. Survival analysis----------------------------------------------------------------------
rm(list = ls())
library(data.table)
library(survminer)
library(survMisc)
library(survival)
library(patchwork)
source("Function.R")

## 2.1 Survival analysis of five bulk data---------------------------------------------------------------------
### load data
Bulk_data <- readRDS("./outdata/Bulk_data.rds")
Surplot_list <- list()

data_demo1 <- data.frame(subtype=Bulk_data$TCGA$predict$Subtype,Bulk_data$TCGA$pdata[,2:7])
data_demo2 <- data.frame(subtype=Bulk_data$ICGC$predict$Subtype,Bulk_data$ICGC$pdata[,2:6])
data_demo3 <- data.frame(subtype=Bulk_data$GSE14520$predict$Subtype,Bulk_data$GSE14520$pdata[,2:8])
data_demo4 <- data.frame(subtype=Bulk_data$GSE116174$predict$Subtype,Bulk_data$GSE116174$pdata[,2:6])
data_demo5 <- data.frame(subtype=Bulk_data$"In-house"$predict$Subtype,Bulk_data$"In-house"$pdata)

##OS
Surplot_list$TCGA <- Surplot(data_demo1,type = "OS",title = "TCGA")
Surplot_list$ICGC <- Surplot(data_demo2,type = "OS",title = "ICGC")
Surplot_list$GSE14520 <- Surplot(data_demo3,type = "OS",title = "GSE14520")
Surplot_list$GSE116174 <- Surplot(data_demo4,type = "OS",title = "GSE116174")
Surplot_list$"In-house" <- Surplot(data_demo5,type = "OS",title = "In-house")

###PFI in TCGA and RFS in GSE14520 and In-house
{
  ### 
  Surplot_list$TCGA_PFI <- Surplot(data_demo1,type = "PFI",title = "TCGA")
  Surplot_list$GSE14520_RFS <- Surplot(data_demo3,type = "RFS",title = "GSE14520")
  Surplot_list$"In-house_RFS" <- Surplot(data_demo5,type = "RFS",title = "In-house")
}
saveRDS(Surplot_list,file = "./outdata/survival/Surplot_list.rds")
### patchwork

pdf("./plot/Surplot_five.pdf",width = 12,height = 12,onefile = F)
  Surplot_list$TCGA$plot + Surplot_list$ICGC$plot + Surplot_list$GSE14520$plot +
    Surplot_list$GSE116174$plot +Surplot_list$`In-house`$plot+
  Surplot_list$TCGA_PFI$plot + Surplot_list$GSE14520_RFS$plot+
    Surplot_list$`In-house_RFS`$plot
dev.off()


## 2.2 Survival analysis of different satges -------------------------------------------------------
rm(list = ls())
library(data.table)
library(survminer)
library(survMisc)
library(survival)
library(patchwork)
source("Function.R")
### load data
Bulk_data <- readRDS("./outdata/Bulk_data.rds")
colnames(Bulk_data$`In-house`$pdata)[10]<-"stage"

###combine data
stage_all <- lapply(Bulk_data, function(f){
  x <- f$pdata
  colnames(x)[1]<-"sample"
  temp <- c("sample","OS","OS.time","stage","PFI","PFI.time","RFS","RFS.time")
  return(data.frame(subtype=f$predict$Subtype,x[,na.omit(match(temp,colnames(x)))]))
})
saveRDS(stage_all,file = "./outdata/survival/stage_all.rds")
stage_all<-readRDS( "./outdata/survival/stage_all.rds")
##plot
satge_plot <- mapply(function(x,y){
  name = c("Stage I-II","Stage III-IV")
  satge_plot<-list()
  for(i in c(1,2)){
      if(i==1){
        temp=c(1,2)
      }else{
        temp=c(3,4)
      }
      data_demo  <- data.frame(x[x$stage %in% temp,])
      for (j in c("OS","RFS","PFI")) {
        if(j %in% colnames(x)){
          satge_plot[paste0(j,"-",name[i])] <- Surplot(data_demo,j,paste0(y," ",name[i]))
        }
      }
  }
  return(satge_plot)
  },stage_all,names(stage_all),SIMPLIFY = F)
##OS all
satge_OS <- do.call(rbind,lapply(stage_all, function(x){
  x[c("subtype","sample","OS","OS.time","stage")]
}))
name = c("Stage I","Stage II","Stage III","Stage IV")
for (i in c(1:4)) {
  data_demo  <- data.frame(satge_OS[satge_OS$stage %in% i,])
  satge_plot$all_os[name[i]] <- Surplot(data_demo,"OS",paste0("OS ",name[i]))
}
###RFS all
satge_RFS <- do.call(rbind,lapply(stage_all, function(x){
  if("RFS" %in% colnames(x)){
      x[c("subtype","sample","RFS","RFS.time","stage")]
  }
}))
name = c("Stage I","Stage II","Stage III","Stage IV")
for (i in c(1:3)) {
  data_demo  <- data.frame(satge_RFS[satge_RFS$stage %in% i,])
  satge_plot$all_rfs[name[i]] <- Surplot(data_demo,"RFS",paste0("RFS ",name[i]))
}
saveRDS(satge_plot,file = "./outdata/survival/satge_plot.rds")

###patchwork
library(patchwork)
pdf("./plot/stage_all_sur.pdf",width = 20,height = 20,onefile = F)
wrap_plots(satge_plot$TCGA$`OS-Stage I-II`,satge_plot$TCGA$`OS-Stage III-IV`,
           satge_plot$ICGC$`OS-Stage I-II`,satge_plot$ICGC$`OS-Stage III-IV`,
           satge_plot$GSE14520$`OS-Stage I-II`,satge_plot$GSE14520$`OS-Stage III-IV`,
           satge_plot$GSE116174$`OS-Stage I-II`,satge_plot$GSE116174$`OS-Stage III-IV`,
           satge_plot$`In-house`$`OS-Stage I-II`,satge_plot$`In-house`$`OS-Stage III-IV`,
           satge_plot$TCGA$`PFI-Stage I-II`,satge_plot$TCGA$`PFI-Stage III-IV`,
           satge_plot$GSE14520$`RFS-Stage I-II`,satge_plot$GSE14520$`RFS-Stage III-IV`,
           satge_plot$`In-house`$`RFS-Stage I-II`,satge_plot$`In-house`$`RFS-Stage III-IV`,
           satge_plot$all_os$`Stage I`,satge_plot$all_os$`Stage II`,satge_plot$all_os$`Stage III`,
           satge_plot$all_rfs$`Stage I`,satge_plot$all_rfs$`Stage II`,satge_plot$all_rfs$`Stage III`
           )
dev.off()

## 2.3 The sample distribution of TNM stage  -----------------------------------------------------------

{
  library(data.table)
  library(ggpubr)
  library(ggstatsplot)
  library(patchwork)
  library(plyr)
  library(ggsci)
  library(ggh4x)
  
  
  satge_all1 <- do.call(rbind,mapply(function(x,y){
    data.frame(dataset=y,x[c("subtype","stage")])
  },stage_all,names(stage_all),SIMPLIFY = F))
  
  a <- data.frame(table(satge_all1$subtype,satge_all1$stage,satge_all1$dataset))
  a$Var4 <-paste0(a$Var1,"-",a$Var3)
  a<- ddply(a,.(Var4),transform,percent=Freq/sum(Freq)*100)
  a$Var2 <- factor(a$Var2,
                   levels = c('1','2','3',"4"),
                   labels = c("I","II","III","IV"))
  a$Var3 <- factor(a$Var3,
                   levels = c('TCGA','ICGC','GSE14520',"GSE116174","In-house"),
                   labels = c("TCGA ****","ICGC ****","GSE14520","GSE116174","In-house ***"))
  p1<-ggplot(a,aes(Var4,percent,fill=Var2))+
    geom_bar(stat="identity",position = position_stack())+
    theme_bw(base_size = 15)+
    theme(panel.grid=element_blank())+
    coord_cartesian(clip = "off")+
    scale_x_discrete(breaks=unique(a$Var4),labels=substring(unique(a$Var4),1,2))+
    scale_y_continuous(labels = scales::percent_format(scale = 1))+
    theme(axis.ticks.x=element_line(color="black"),
          axis.ticks.y =element_line(color="black"),
          axis.line.y=element_line(color="black"),
          axis.line.x=element_line(color="black"))+
    font("xy.text", size = 14,color = "black")+
    font("ylab", size = 15,color = "black")+
    labs(y="Percent%",color = "Stage",fill= "Stage")+
    font("legend.text", size = 14,face = "plain")+
    font("legend.title", size = 14)+
    scale_fill_nejm()+
    theme(axis.title.x =element_blank())+
    facet_nested(~ Var3,scales = "free")+
    theme(strip.text.x = element_text(size = 14),
          strip.placement = "outside")+
    theme(strip.background.x = element_rect(fill = "grey85", colour = "grey85"),
          panel.border=element_rect(colour ="black"),
          axis.text.x  = element_text(colour = c("#CD1818","#0F2C67")))
  
  pdf("./plot/stage_bar.pdf",height = 4,width = 8)
  p1
  dev.off()
}

## 2.4 Univariate and multivariate Cox regression analyses -----------------------------------------------------------------
library(forestplot)
library(survival)
Bulk_data <- readRDS("E:/project/HCC_Immune_Subtype/outdata/Bulk_data.rds")
colnames(Bulk_data$`In-house`$pdata)[c(4,10)]<-c("age2","stage")

cox_result <- lapply(Bulk_data, function(x){
  cox_s_m(x$pdata,x$predict$Subtype)
})

save(cox_result,file ="./outdata/survival/cox_result.Rdata")
cox_result_com <- lapply(cox_result, function(x){
  single <- x$single[,-2]
  colnames(single)<-c("HR","lower95","upper95","Pvalue","U-Hazard Ratio(95%CI)")
  single$Pvalue <- ifelse(single$Pvalue<0.001,"p<0.001",single$Pvalue)
  multi <- x$multi[,-2]
  colnames(multi)<-c("HR","lower95","upper95","Pvalue","U-Hazard Ratio(95%CI)")
  multi$Pvalue <- ifelse(multi$Pvalue<0.001,"p<0.001",multi$Pvalue)
  cbind(single,multi)
})
library(openxlsx)
write.xlsx(cox_result_com,file = "./outdata/survival/cox_result_com.xlsx")

#####画图
TCGA<-read.xlsx("./outdata/survival/cox_result.xlsx",sheet =1)
A1<-rbind(colnames(TCGA),TCGA)
p1<-forestplot(labeltext = as.matrix(A1[,c(1,6,5,11,10)]),
               #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
               mean = cbind(as.numeric(A1[,2]),as.numeric(A1[,7])), #设置均值
               lower = cbind(as.numeric(A1[,3]),as.numeric(A1[,8])), #设置均值的lowlimits限
               upper = cbind(as.numeric(A1[,4]),as.numeric(A1[,9])), #设置均值的uplimits限
               zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
               boxsize = 0.1, #设置点估计的方形大小
               lineheight = unit(15,'mm'),#设置图形中的行距
               colgap = unit(4.5,'mm'),#设置图形中的列间距
               lwd.zero = 2,#设置参考线的粗细
               lwd.ci = 2,#设置区间估计线的粗细
               col=fpColors(box=c('#CD1818','#0F2C67'),lines = c('black','black')),
               lwd.xaxis=2,#设置X轴线的粗细
               lty.ci = "solid",
               graph.pos = 4,
               xticks=c(0,1,4,8,12),
               is.summary=c(T,F,F,F,F,F),
               graphwidth = unit(.15,"npc"),
               hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                               "2" = gpar(lty=2),
                               "7" = gpar(lty=1,lwd=2)),
               legend = c("Univariate model", "Multivariable model"),
               txt_gp=fpTxtGp(
                 label=gpar(cex=1),
                 ticks=gpar(cex=1), 
                 xlab=gpar(cex=1.5), 
                 title=gpar(cex=2)))

ICGC<-read.xlsx("./outdata/survival/cox_result.xlsx",sheet =2)
A1<-rbind(colnames(ICGC),ICGC)
p2<-forestplot(labeltext = as.matrix(A1[,c(1,6,5,11,10)]),
               #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
               mean = cbind(as.numeric(A1[,2]),as.numeric(A1[,7])), #设置均值
               lower = cbind(as.numeric(A1[,3]),as.numeric(A1[,8])), #设置均值的lowlimits限
               upper = cbind(as.numeric(A1[,4]),as.numeric(A1[,9])), #设置均值的uplimits限
               zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
               boxsize = 0.1, #设置点估计的方形大小
               lineheight = unit(15,'mm'),#设置图形中的行距
               colgap = unit(4.5,'mm'),#设置图形中的列间距
               lwd.zero = 2,#设置参考线的粗细
               lwd.ci = 2,#设置区间估计线的粗细
               col=fpColors(box=c('#CD1818','#0F2C67'),lines = c('black','black')),
               lwd.xaxis=2,#设置X轴线的粗细
               lty.ci = "solid",
               graph.pos = 4,
               xticks=c(0,1,4,8),
               is.summary=c(T,F,F,F,F,F),
               graphwidth = unit(.15,"npc"),
               hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                               "2" = gpar(lty=2),
                               "6" = gpar(lty=1,lwd=2)),
               legend = c("Univariate model", "Multivariable model"),
               txt_gp=fpTxtGp(
                 label=gpar(cex=1),
                 ticks=gpar(cex=1), 
                 xlab=gpar(cex=1.5), 
                 title=gpar(cex=2)))

GSE14520<-read.xlsx("./outdata/survival/cox_result.xlsx",sheet =3)
A1<-rbind(colnames(GSE14520),GSE14520)
p3<-forestplot(labeltext = as.matrix(A1[,c(1,6,5,11,10)]),
               #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
               mean = cbind(as.numeric(A1[,2]),as.numeric(A1[,7])), #设置均值
               lower = cbind(as.numeric(A1[,3]),as.numeric(A1[,8])), #设置均值的lowlimits限
               upper = cbind(as.numeric(A1[,4]),as.numeric(A1[,9])), #设置均值的uplimits限
               zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
               boxsize = 0.1, #设置点估计的方形大小
               lineheight = unit(15,'mm'),#设置图形中的行距
               colgap = unit(4.5,'mm'),#设置图形中的列间距
               lwd.zero = 2,#设置参考线的粗细
               lwd.ci = 2,#设置区间估计线的粗细
               col=fpColors(box=c('#CD1818','#0F2C67'),lines = c('black','black')),
               lwd.xaxis=2,#设置X轴线的粗细
               lty.ci = "solid",
               graph.pos = 4,
               xticks=c(0,1,2,3,4),
               is.summary=c(T,F,F,F,F,F),
               graphwidth = unit(.15,"npc"),
               hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                               "2" = gpar(lty=2),
                               "6" = gpar(lty=1,lwd=2)),
               legend = c("Univariate model", "Multivariable model"),
               txt_gp=fpTxtGp(
                 label=gpar(cex=1),
                 ticks=gpar(cex=1), 
                 xlab=gpar(cex=1.5), 
                 title=gpar(cex=2)))

GSE116174<-read.xlsx("./outdata/survival/cox_result.xlsx",sheet =4)
A1<-rbind(colnames(GSE116174),GSE116174)
p4<-forestplot(labeltext = as.matrix(A1[,c(1,6,5,11,10)]),
               #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
               mean = cbind(as.numeric(A1[,2]),as.numeric(A1[,7])), #设置均值
               lower = cbind(as.numeric(A1[,3]),as.numeric(A1[,8])), #设置均值的lowlimits限
               upper = cbind(as.numeric(A1[,4]),as.numeric(A1[,9])), #设置均值的uplimits限
               zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
               boxsize = 0.1, #设置点估计的方形大小
               lineheight = unit(15,'mm'),#设置图形中的行距
               colgap = unit(4.5,'mm'),#设置图形中的列间距
               lwd.zero = 2,#设置参考线的粗细
               lwd.ci = 2,#设置区间估计线的粗细
               col=fpColors(box=c('#CD1818','#0F2C67'),lines = c('black','black')),
               lwd.xaxis=2,#设置X轴线的粗细
               lty.ci = "solid",
               graph.pos = 4,
               xticks=c(0,1,4,8),
               is.summary=c(T,F,F,F,F,F),
               graphwidth = unit(.15,"npc"),
               hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                               "2" = gpar(lty=2),
                               "6" = gpar(lty=1,lwd=2)),
               legend = c("Univariate model", "Multivariable model"),
               txt_gp=fpTxtGp(
                 label=gpar(cex=1),
                 ticks=gpar(cex=1), 
                 xlab=gpar(cex=1.5), 
                 title=gpar(cex=2)))

Inhouse<-read.xlsx("./outdata/survival/cox_result.xlsx",sheet =5)
A1<-rbind(colnames(Inhouse),Inhouse)
p5<-forestplot(labeltext = as.matrix(A1[,c(1,6,5,11,10)]),
               #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
               mean = cbind(as.numeric(A1[,2]),as.numeric(A1[,7])), #设置均值
               lower = cbind(as.numeric(A1[,3]),as.numeric(A1[,8])), #设置均值的lowlimits限
               upper = cbind(as.numeric(A1[,4]),as.numeric(A1[,9])), #设置均值的uplimits限
               zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
               boxsize = 0.1, #设置点估计的方形大小
               lineheight = unit(15,'mm'),#设置图形中的行距
               colgap = unit(4.5,'mm'),#设置图形中的列间距
               lwd.zero = 2,#设置参考线的粗细
               lwd.ci = 2,#设置区间估计线的粗细
               col=fpColors(box=c('#CD1818','#0F2C67'),lines = c('black','black')),
               lwd.xaxis=2,#设置X轴线的粗细
               lty.ci = "solid",
               graph.pos = 4,
               xticks=c(0,1,4,8),
               is.summary=c(T,F,F,F,F,F),
               graphwidth = unit(.15,"npc"),
               hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                               "2" = gpar(lty=2),
                               "6" = gpar(lty=1,lwd=2)),
               legend = c("Univariate model", "Multivariable model"),
               txt_gp=fpTxtGp(
                 label=gpar(cex=1),
                 ticks=gpar(cex=1), 
                 xlab=gpar(cex=1.5), 
                 title=gpar(cex=2)))

####所有29对标志合在一起
genepair<-read.xlsx("./outdata/survival/cox_result.xlsx",sheet =6)
A1<-rbind(colnames(genepair),genepair)
p6<-forestplot(labeltext = as.matrix(A1[,c(1,6,5,11,10)]),
               #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
               mean = cbind(as.numeric(A1[,2]),as.numeric(A1[,7])), #设置均值
               lower = cbind(as.numeric(A1[,3]),as.numeric(A1[,8])), #设置均值的lowlimits限
               upper = cbind(as.numeric(A1[,4]),as.numeric(A1[,9])), #设置均值的uplimits限
               zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
               boxsize = 0.1, #设置点估计的方形大小
               lineheight = unit(15,'mm'),#设置图形中的行距
               colgap = unit(4.5,'mm'),#设置图形中的列间距
               lwd.zero = 2,#设置参考线的粗细
               lwd.ci = 2,#设置区间估计线的粗细
               col=fpColors(box=c('#CD1818','#0F2C67'),lines = c('black','black')),
               lwd.xaxis=2,#设置X轴线的粗细
               lty.ci = "solid",
               graph.pos = 4,
               xticks=c(0,1,4,8,12),
               is.summary=c(T,F,F,F,F,F),
               graphwidth = unit(.15,"npc"),
               hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                               "2" = gpar(lty=2),
                               "7" = gpar(lty=1,lwd=2)),
               legend = c("Univariate model", "Multivariable model"),
               txt_gp=fpTxtGp(
                 label=gpar(cex=1),
                 ticks=gpar(cex=1), 
                 xlab=gpar(cex=1.5), 
                 title=gpar(cex=2)))
###
pdf("./plot/cox/cox_TCGA.pdf",width=15,height=6)
p1
dev.off()
pdf("./plot/cox/cox_ICGC.pdf",width=15,height=6)
p2
dev.off()
pdf("./plot/cox/cox_GSE14520.pdf",width=15,height=6)
p3
dev.off()
pdf("./plot/cox/cox_GSE116174.pdf",width=15,height=6)
p4
dev.off()
pdf("./plot/cox/cox_Inhouse.pdf",width=15,height=6)
p5
dev.off()
pdf("./plot/cox/cox_pair.pdf",width=15,height=6)
p6
dev.off()

## 2.5 ROC -----------------------------------------------------------------
### roc
ROC_hiplot <- list()
name <- c("TCGA","ICGC","GSE14520","GSE116174","In-house")
for (i in 1:5) {
  ROC_hiplot[[name[i]]]<-data.frame(time=Bulk_data[[name[i]]]$pdata$OS.time,
                                    cens=Bulk_data[[name[i]]]$pdata$OS,
                                    ratio=Bulk_data[[name[i]]]$predict$Ratio)
  colnames(ROC_hiplot[[name[i]]])[3]<-"ratio"
}
saveRDS(ROC_hiplot,file = "./outdata/survival/dataROC.rds")

library(survivalROC)
dataset <- c("TCGA","ICGC","GSE14520","GSE116174","In-house")
roc_plot <- list()
for (j in dataset) {
  TCGA <- ROC_hiplot[[j]]
  p_time <- c(1,3,5)*365#设置相应的时间点
  time_roc_res <- list(FP=NULL,TP=NULL,AUC=NULL)
  for (i in p_time){
    ROC.alb<-survivalROC(Stime=TCGA$time,
                         status=TCGA$cens,
                         marker=TCGA$ratio,
                         predict.time = i,
                         method="KM")
    time_roc_res$FP <- cbind(time_roc_res$FP,ROC.alb$FP)
    time_roc_res$TP <- cbind(time_roc_res$TP,ROC.alb$TP)
    time_roc_res$AUC <- cbind(time_roc_res$AUC,ROC.alb$AUC)
  }
  time_ROC_df <- data.frame(
    TP_3year = time_roc_res$TP[,1],
    FP_3year = time_roc_res$FP[,1],
    TP_5year = time_roc_res$TP[,2],
    FP_5year = time_roc_res$FP[,2],
    TP_10year = time_roc_res$TP[,3],
    FP_10year = time_roc_res$FP[,3]
  )
  p1 <- ggplot(data = time_ROC_df) +
    geom_line(aes(x = FP_3year, y = TP_3year), size = 1.2, color = "#BC3C29FF") +
    geom_line(aes(x = FP_5year, y = TP_5year), size = 1.2, color = "#0072B5FF") +
    geom_line(aes(x = FP_10year, y = TP_10year), size = 1.2, color = "#E18727FF") +
    geom_abline(slope = 1, intercept = 0, color = "grey75", size = 1, linetype = 2) +
    theme_bw() +
    annotate("text",
             x = 0.7, y = 0.25, size = 5,
             label = paste0("AUC at 1 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
    ) +
    annotate("text",
             x = 0.7, y = 0.15, size = 5,
             label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
    ) +
    annotate("text",
             x = 0.7, y = 0.05, size = 5,
             label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
    ) +
    labs(x = "False positive rate", y = "True positive rate") +
    theme(
      axis.text = element_text(size = 14, color = "black"),
      axis.title.x = element_text(size = 15, color = "black", margin = margin(c(15, 0, 0, 0))),
      axis.title.y = element_text(size = 15, color = "black", margin = margin(c(0, 15, 0, 0)))
    )+
    ggtitle(j)
  roc_plot[[j]] <- p1
}

library(patchwork)
pdf("./plot/ROC.pdf",width=12,height=8)
wrap_plots(roc_plot)
dev.off()


# 3. The immune landscape----------------------------------------------------------------------
## 3.1 Immune infiltration -------------------------------------------------
rm(list = ls())
library(data.table)
library(sqldf)
library(dplyr)
library(survival)
library(foreach)
library(utils)
library(survminer)
library(survMisc)
library(tibble)
library(vegan)
library(patchwork)
#
savedir<-"./outdata/immune_microenvironment/"
source("Function.R")
##
library(immunedeconv)
Bulk_exp <- readRDS("./data/Bulk_exp.rds")
Bulk_data <- readRDS("./outdata/Bulk_data.rds")
##
TCGA<-as.data.frame(egg_us(Bulk_exp$TCGA))
TCGA_inf<-Im.inf(TCGA,array=FALSE)
rm(TCGA)
###
ICGC<-as.data.frame(egg_us(Bulk_exp$ICGC))
ICGC_inf<-Im.inf(ICGC,array=FALSE)
rm(ICGC)
##
GSE14520<-as.data.frame(egg_us(Bulk_exp$GSE14520))
GSE14520_inf<-Im.inf(GSE14520,array=TRUE)
rm(GSE14520)
##
GSE116174<-as.data.frame(egg_us(Bulk_exp$GSE116174))
GSE116174_inf<-Im.inf(GSE116174,array=TRUE)
rm(GSE116174)
##
In_house<-as.data.frame(egg_us(Bulk_exp$"In-house"))
In_house_inf<-Im.inf(In_house,array=F)
rm(In_house)

inf<-list(TCGA=TCGA_inf,ICGC=ICGC_inf,GSE14520=GSE14520_inf,GSE116174=GSE116174_inf,
          "In-house"=In_house_inf)
save(inf,file = paste0(savedir,"inf.Rdata"))
load(paste0(savedir,"inf.Rdata"))
temp<-lapply(inf,function(x){
  z<-as.list(names(x))
  mapply(y=x,function(y,z){
    rownames(y)<-paste(rownames(y),sep="_",z)
    return(list(t(y)))
  },z)
})

label_o<-list(TCGA=Bulk_data$TCGA$predict$Subtype,
              ICGC=Bulk_data$ICGC$predict$Subtype,
              GSE14520=Bulk_data$GSE14520$predict$Subtype,
              GSE116174=Bulk_data$GSE116174$predict$Subtype,
              "In-house" = Bulk_data$"In-house"$predict$Subtype)
inf.cbind<-lapply(temp,function(x){as.data.frame(do.call('cbind',lapply(x,function(x) x)))})
inf.cbind<-mapply(function(x,y){data.frame(x,y)},y=inf.cbind,x=label_o)
save(inf.cbind,file = paste0(savedir,"inf.cbind.Rdata"))
load(paste0(savedir,"inf.cbind.Rdata"))

##cibersort
source("./code/CIBERSORT.R")
##
ciber<-lapply(Bulk_exp, function(x){
  print("1")
  if(length(grep("GSM",x))==0){
    TME.results = CIBERSORT("./data/LM22.txt", 
                            x , 
                            perm = 1000, 
                            QN = T)
  }
  else{
    TME.results = CIBERSORT("./data/LM22.txt", 
                            x , 
                            perm = 1000, 
                            QN = F)
  }})

ciber<-lapply(ciber,function(y){
  colnames(y)<-paste(colnames(y),sep="_","ciber")
  y<-as.data.frame(y)
  return(y)
})
save(ciber,file = paste0(savedir,"ciber.Rdata"))
load(paste0(savedir,"ciber.Rdata"))
##
inf.five<-mapply(y=inf.cbind,function(y,z){
  z<-z[,1:22]
  a<-data.frame(y,z)
  return(a)
},z=ciber)
save(inf.five,file = paste0(savedir,"inf.five.Rdata"))
load(paste0(savedir,"inf.five.Rdata"))

######## 2.SSGSEA
##SSGSEA
library(GSVA)
set_30594216<-freadlist("./data/30594216set.txt")
set_24<-freadlist("./data/APM_24IMMUNE.txt")
set_30594216$`CD8+ T cells`<-NULL
set_30594216$`T helper cells`<-NULL
set_24$APM<-NULL
set_24$Treg<-NULL
set_24$pDC<-NULL

imm.29 <- lapply(Bulk_exp, function(x){
  x <- egg_us(x)
  out <- gsva(x,set_30594216,method='ssgsea')
}) 
##
imm.29<-lapply(imm.29,function(y){
  y<-t(y)
  y<-as.data.frame(y)
  colnames(y)<-paste(colnames(y),sep="_","29")
  return(y)
})
##
imm.24 <- lapply(Bulk_exp, function(x){
  x <- egg_us(x)
  out <- gsva(x,set_24,method='ssgsea')
}) 
##
imm.24<-lapply(imm.24,function(y){
  y<-t(y)
  y<-as.data.frame(y)
  colnames(y)<-paste(colnames(y),sep="_","24")
  return(y)
})
###
inf.ssGSEA<-mapply(function(x,y){data.frame(x,y)},x=imm.29,y=imm.24,SIMPLIFY = FALSE)
save(inf.ssGSEA,file =paste0(savedir,"inf.ssGSEA.Rdata"))
load(paste0(savedir,"inf.ssGSEA.Rdata"))
##########3.TIDE
TIDE<-list()
TIDE$TCGA <- fread("./data/TCGA_TIDE.csv",data.table = F)[,c(1,13:15)]
TIDE$TCGA <- TIDE$TCGA[match(substring(rownames(ciber$TCGA),1,16),substring(TIDE$TCGA[,1],1,16)),]
TIDE$TCGA <- egg_us(TIDE$TCGA)
colnames(TIDE$TCGA)<-c("MDSC_TIDE","CAF_TIDE","M2TAM_TIDE")

TIDE$ICGC <- fread("./data/ICGC_TIDE.csv",data.table = F)[,c(1,13:15)]
TIDE$ICGC <- TIDE$ICGC[match(rownames(ciber$ICGC),TIDE$ICGC[,1]),]
TIDE$ICGC <- egg_us(TIDE$ICGC)
colnames(TIDE$ICGC)<-c("MDSC_TIDE","CAF_TIDE","M2TAM_TIDE")

TIDE$GSE14520 <- fread("./data/GSE14520_TIDE.csv",data.table = F)[,c(1,13:15)]
TIDE$GSE14520 <- TIDE$GSE14520[match(rownames(ciber$GSE14520),TIDE$GSE14520[,1]),]
TIDE$GSE14520 <- egg_us(TIDE$GSE14520)
colnames(TIDE$GSE14520)<-c("MDSC_TIDE","CAF_TIDE","M2TAM_TIDE")

TIDE$GSE116174 <- fread("./data/GSE116174_TIDE.csv",data.table = F)[,c(1,13:15)]
TIDE$GSE116174 <- TIDE$GSE116174[match(rownames(ciber$GSE116174),TIDE$GSE116174[,1]),]
TIDE$GSE116174 <- egg_us(TIDE$GSE116174)
colnames(TIDE$GSE116174)<-c("MDSC_TIDE","CAF_TIDE","M2TAM_TIDE")

TIDE$"In-house" <- fread("./data/Inhouse_TIDE.csv",data.table = F)[,c(1,13:15)]
TIDE$"In-house" <- TIDE$"In-house"[match(rownames(ciber$"In-house"),TIDE$"In-house"[,1]),]
TIDE$"In-house" <- egg_us(TIDE$"In-house")
colnames(TIDE$"In-house")<-c("MDSC_TIDE","CAF_TIDE","M2TAM_TIDE")
##
###CIBERSORT
Lymphocytes_C<-c("B.cells.naive_ciber","B.cells.memory_ciber","T.cells.CD4.naive_ciber",
                 "T.cells.CD4.memory.resting_ciber","T.cells.CD4.memory.activated_ciber",
                 "T.cells.follicular.helper_ciber","T.cells.regulatory..Tregs._ciber",
                 "T.cells.gamma.delta_ciber","T.cells.CD8_ciber","NK.cells.resting_ciber",
                 "NK.cells.activated_ciber","Plasma.cells_ciber")
Macrophages_C<-c("Monocytes_ciber","Macrophages.M0_ciber",
                 "Macrophages.M1_ciber","Macrophages.M2_ciber")
Macrophages_Q<-c("Monocyte_quantiseq","Macrophage.M1_quantiseq","Macrophage.M2_quantiseq")
##
CQ_co<-list(Lymphocytes_C=Lymphocytes_C,Macrophages_C=Macrophages_C,
            Macrophages_Q=Macrophages_Q)
all_cq<-lapply(inf.five, function(y){
  test<-as.data.frame(sapply(CQ_co, function(x){
    apply(y[,x,drop=FALSE],1,sum)
  }))
  test<-data.frame(subtype=y[,1,drop=FALSE],test)
})

##########合并替换
inf.all<-mapply(function(x,y,z,q){data.frame(x,y,z,q)},x=inf.five,y=inf.ssGSEA,
                z=TIDE,q=all_cq,SIMPLIFY = FALSE)
inf.all<-lapply(inf.all, function(x){
  colnames(x)[1]<-"subtype"
  return(x)
})
save(inf.all,file=paste0(savedir,"inf.all.Rdata"))
library(readxl)
data<-xlsx.list("./data/inf.xlsx")
TCGA<-imm.select(inf.all,"TCGA",data)
ICGC<-imm.select(inf.all,"ICGC",data)
GSE14520<-imm.select(inf.all,"GSE14520",data)
GSE116174<-imm.select(inf.all,"GSE116174",data)
Inhouse <- imm.select(inf.all,"In-house",data)
inf.summary<-list(TCGA=TCGA,ICGC=ICGC,GSE14520=GSE14520,GSE116174=GSE116174,"In-house"=Inhouse)
save(inf.summary,file =paste0(savedir,"inf.summary.Rdata"))
load(paste0(savedir,"inf.summary.Rdata"))
rm(TCGA,ICGC,GSE14520,GSE116174,Inhouse)


##
data<-lapply(data,function(x){
  colnames(x)[2]<-"variable"
  return(x)})
inf.name<-do.call("rbind",lapply(data,function(x){x[,2,drop=FALSE]}))
save(inf.name,file =paste0(savedir,"inf.name"))
rm(data)
####
inf.summary1<-lapply(inf.summary, function(x){x$total})
temp<-inf.name
for (i in 1:5) {
  temp<-full_join(temp,inf.summary1[[i]][,c(1,8)],by="variable")
}
temp1<-which(apply(temp,1,function(x){
  sum(is.na(x))
})==5)
temp<-temp[-temp1,]
colnames(temp)<-c("cell","TCGA","ICGC","GSE14520","GSE116174","In-house")
temp<-as.data.frame(temp)
library(tidyverse)

####1:55,56:98,101:102
data<-reshape2::melt(temp[c(55:1),])
data<-data.frame(data,color=((data$value>0)-0.5)*2)
data<-replace_egg(data,4,from1 = -1,to1 = "C2",from2 = 1,to2 = "C1")
data1<-data %>% mutate(pvalue = case_when(
  abs(data$value) <0.05 ~ "p<0.05",
  abs(data$value) >=0.05 ~"p>0.05",
  is.na(data$value)==TRUE  ~  "na" )
)
data1$pvalue<-factor(data1$pvalue,levels = c("p>0.05","p<0.05","na"))
data1$color<-factor(data1$color,levels = c("C1","C2"))
data1$cell<-factor(data1$cell,levels = as.character(unique(data1$cell)))
data1<-na.omit(data1)
p1<-ggplot(data1, aes(x = variable, y = cell, size=pvalue,color=color)) + 
  geom_point()+
  cowplot::theme_cowplot() +
  theme(axis.line = element_blank()) +
  #theme(panel.border = element_rect(fill=NA,color="black"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=1)) +
  ylab(NULL) +
  theme(axis.ticks = element_blank(),axis.title.x=element_blank()) +
  scale_y_discrete(position = "left")+
  scale_color_manual(values = c("C1"="#CD1818","C2"="#0F2C67"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  theme(legend.text = element_text(size = 14))
p1

data<-reshape2::melt(temp[c(100:56),])
data<-data.frame(data,color=((data$value>0)-0.5)*2)
data<-replace_egg(data,4,from1 = -1,to1 = "C2",from2 = 1,to2 = "C1")
data1<-data %>% mutate(pvalue = case_when(
  abs(data$value) <0.05 ~ "p<0.05",
  abs(data$value) >=0.05 ~"p>0.05",
  is.na(data$value)==TRUE  ~  "na" )
)
data1$pvalue<-factor(data1$pvalue,levels = c("p>0.05","p<0.05","na"))
data1$color<-factor(data1$color,levels = c("C1","C2"))
data1$cell<-factor(data1$cell,levels = as.character(unique(data1$cell)))
data1<-na.omit(data1)
p2<-ggplot(data1, aes(x = variable, y = cell, size=pvalue,color=color)) + 
  geom_point()+
  cowplot::theme_cowplot() +
  theme(axis.line = element_blank()) +
  #theme(panel.border = element_rect(fill=NA,color="black"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=1)) +
  ylab(NULL) +
  theme(axis.ticks = element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank()) +
  scale_y_discrete(position = "right")+
  scale_color_manual(values = c("C1"="#CD1818","C2"="#0F2C67"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  theme(legend.text = element_text(size = 14))
p2

data<-reshape2::melt(temp[c(113:101),])
data<-data.frame(data,color=((data$value>0)-0.5)*2)
data<-replace_egg(data,4,from1 = -1,to1 = "C2",from2 = 1,to2 = "C1")
data1<-data %>% mutate(pvalue = case_when(
  abs(data$value) <0.05 ~ "p<0.05",
  abs(data$value) >=0.05 ~"p>0.05",
  is.na(data$value)==TRUE  ~  "na" )
)
data1$pvalue<-factor(data1$pvalue,levels = c("p>0.05","p<0.05","na"))
data1$color<-factor(data1$color,levels = c("C1","C2"))
data1$cell<-factor(data1$cell,levels = as.character(unique(data1$cell)))
data1<-na.omit(data1)
p3<-ggplot(data1, aes(x = variable, y = cell, size=pvalue,color=color)) + 
  geom_point()+
  cowplot::theme_cowplot() +
  theme(axis.line = element_blank()) +
  #theme(panel.border = element_rect(fill=NA))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=1)) +
  ylab(NULL) +
  theme(axis.ticks = element_blank(),axis.title.x=element_blank()) +
  scale_y_discrete(position = "right")+
  scale_color_manual(values = c("C1"="#CD1818","C2"="#0F2C67"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  theme(legend.text = element_text(size = 14))
p3


library(patchwork)
p4<-p1+(p2/p3 + plot_layout(heights = c(45,13)))+plot_layout(guides = 'collect')

pdf("./plot/inf.pdf",width=8,height=14)
p4
dev.off()

## 3.2 10 gene sets -----------------------------------------------------------
savedir<-"./outdata/immune_microenvironment/"
library(GSVA)
library(data.table)
library(ggplot2)
library(ggpubr)
source("Function.R")
geneset <- freadlist("./data/Geneset.txt")
Bulk_exp<-readRDS("./data/Bulk_exp.rds")
predict_data <- readRDS("./outdata/predict_data.rds")
out <- list() 
name <- c("TCGA","ICGC","GSE14520","GSE116174","In-house")
for (i in 1:5) {
  data <- Bulk_exp[[name[i]]]
  data <- egg_us(data)
  eggtest <- gsva(data,geneset,method='ssgsea')
  temp <- data.frame(predict_data[[name[i]]]$Subtype,t(eggtest))  
  temp1 <- reshape2::melt(temp)
  temp1 <- data.frame(temp1,dataset=name[i])
  out[[name[i]]] <- temp1
}
Imm_feature <- do.call(rbind,lapply(out,function(x) x))
colnames(Imm_feature)[1]<-"subtype"
Imm_feature$dataset <- factor(Imm_feature$dataset,levels = c("TCGA","ICGC","GSE14520","GSE116174","In-house"))
Imm_feature$subtype <- factor(Imm_feature$subtype,levels = c("C1","C2"))
saveRDS(Imm_feature,file =paste0(savedir,"Imm_feature.rds"))
name <- as.character(unique(Imm_feature$variable))
name1 <- names(geneset)
plot_imm <- list() 
for (i in 1:10) {
  test<-Imm_feature[Imm_feature$variable == name[i],]
  test$variable <- as.character(test$variable)
  plot_imm[[name1[i]]] <- plot.violin(test,x=dataset,y=value,
                                      fill = subtype,
                                      palette.fill=c("#CD1818","#0F2C67"),
                                      palette.color = c("#CD1818","#0F2C67"),
                                      color = "white",
                                      fill.b = "white",
                                      color.b = subtype,
                                      group=subtype,
                                      alpha.v = 0.7,
                                      notch.b=T,
                                      width.b=0.5,
                                      add = c("boxplot","violinplot","shadow"),#,"shadow"
                                      shadow.col = c("grey40","white"),
                                      x_label="",
                                      y_label = "Score",
                                      theme = theme_bw(base_size = 16))+
    theme(axis.text.x = element_blank(),   ## 删去所有刻度标签
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 1,unit = 'cm'))+
    theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5))+
    coord_cartesian(clip = 'off')
}
plot_imm1 <- lapply(plot_imm, function(x){
  x+facet_wrap(~variable,scales = "free_x")+
    theme(strip.text.x = element_text(size = 15, colour = "black")) + # 设置分面的字字体大小、颜色、背景、边框，
    theme(strip.background.x = element_rect(colour = NA))+
    geom_vline(xintercept = c(1.5,2.5,3.5,4.5), lwd = 1, colour = "grey70")
})
library(jjAnno)
library(dittoSeq)
dittoColors()[1:5]
plot_CTA <- annoRect(object = plot_imm$CTA,
                     annoPos = 'botomn',
                     xPosition = c(1:5),
                     pFill=dittoColors()[1:5],
                     pCol=dittoColors()[1:5],
                     yPosition = c(-0.46,-0.5),
                     rectWidth = 1)
plot_CTA <- plot_CTA+facet_wrap(~variable,scales = "free_x")+
  theme(strip.text.x = element_text(size = 15, colour = "black")) + # 设置分面的字字体大小、颜色、背景、边框，
  theme(strip.background.x = element_rect(colour = NA))+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), lwd = 1, colour = "grey70")

plot_IFN2 <- annoRect(object = plot_imm$`Type II IFN Reponse`,
                      annoPos = 'botomn',
                      xPosition = c(1:5),
                      pFill=dittoColors()[1:5],
                      pCol=dittoColors()[1:5],
                      yPosition = c(-0.25,-0.4),
                      rectWidth = 1)
plot_IFN2 <- plot_IFN2+facet_wrap(~variable,scales = "free_x")+
  theme(strip.text.x = element_text(size = 15, colour = "black")) + # 设置分面的字字体大小、颜色、背景、边框，
  theme(strip.background.x = element_rect(colour = NA))+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), lwd = 1, colour = "grey70")
plot_IFN2


library(patchwork)
pdf(file = "./plot/left_imm.pdf",width =5.5 ,height = 12.5)
plot_imm1$Glycloysis /
  plot_imm1$Hypoxia /
  plot_imm1$Proliferation /
  plot_imm1$Stem /
  plot_CTA
dev.off()

pdf(file = "./plot/right_imm.pdf",width =5.5 ,height = 12.5)
(plot_imm1$Immunosuppression+scale_y_continuous(position = "right")) /
  (plot_imm1$Cytotoxicity+scale_y_continuous(position = "right")) /
  (plot_imm1$`Pro-Inflammation`+scale_y_continuous(position = "right")) /
  (plot_imm1$`Type I IFN Reponse`+scale_y_continuous(position = "right")) /
  (plot_IFN2+scale_y_continuous(position = "right"))
dev.off()


## 3.3 Network -------------------------------------------------------------
###拼表达谱
load("E:/project/HCC_Immune_Subtype/outdata/network/data/node_all.Rdata")
load("E:/project/HCC_Immune_Subtype/outdata/immune_microenvironment/ciber.Rdata")
load("E:/project/HCC_Immune_Subtype/outdata/network/data/all_int.Rdata")
Bulk_exp <- readRDS("E:/project/HCC_Immune_Subtype/data/Bulk_exp.rds")
ciber_cbind<-list(CD19..B.cells=c("B cells naive_ciber","B cells memory_ciber"),
                  CD4..T.cells =c("T cells CD4 naive_ciber","T cells CD4 memory resting_ciber","T cells CD4 memory activated_ciber"),
                  CD8..T.cells=c("T cells CD8_ciber"),
                  Eosinophils=c("Eosinophils_ciber"),
                  Macrophage.Monocyte.derived=c("Macrophages M0_ciber","Macrophages M1_ciber","Macrophages M2_ciber"),
                  Mast.cells=c("Mast cells resting_ciber","Mast cells activated_ciber"),
                  NK.cells=c("NK cells resting_ciber","NK cells activated_ciber"),
                  Neutrophils=c("Neutrophils_ciber"),
                  Dendritic_cells=c("Dendritic cells resting_ciber","Dendritic cells activated_ciber"))
###
network <- list()
network$TCGA <- network_ratio("TCGA",Bulk_exp$TCGA,predict_data$TCGA,all_int,ciber,ciber_cbind)
network$ICGC <- network_ratio("ICGC",Bulk_exp$ICGC,predict_data$ICGC,all_int,ciber,ciber_cbind)
network$GSE14520 <- network_ratio("GSE14520",Bulk_exp$GSE14520,predict_data$GSE14520,all_int,ciber,ciber_cbind)
network$GSE116174 <- network_ratio("GSE116174",Bulk_exp$GSE116174,predict_data$GSE116174,all_int,ciber,ciber_cbind)
network$"In-house" <- network_ratio("In-house",Bulk_exp$"In-house",predict_data$"In-house",all_int,ciber,ciber_cbind)
#保存
save(network,file = "./outdata/network/network.Rdata")
load("./outdata/network/network.Rdata")
test<-lapply(network,function(x){
  p=0.667
  p1=2
  high<-x$high[((x$high$From.ratio > p) & (x$high$To.ratio > p) & (x$high$concordance > p1)),]
  low<-x$low[((x$low$From.ratio > p) & (x$low$To.ratio > p) & (x$low$concordance > p1)),]
  out<-list(high=high,low=low)
  return(out)
}
)

library(tidyverse)
library(plyr)
high<-lapply(test,function(x){x$high})
low<-lapply(test,function(x){x$low})
all_filter<-list(high=high,low=low)

all_filter_three<-lapply(all_filter, function(x){
  x<-rbind(x$TCGA,x$ICGC,x$"In-house")
  x<-x[,1:2]
  out<-ddply(x,.(From,To),nrow)
  return(out)
})

all_filter1 <- all_filter_three
node_high <- all_filter1$high[all_filter1$high$V1>=2,-3]
node_low <- all_filter1$low[all_filter1$low$V1>=2,-3]
high_label<-node_all[match(unique(c(node_high[,1],node_high[,2])),node_all$Node),]
low_label<-node_all[match(unique(c(node_low[,1],node_low[,2])),node_all$Node),]


write.csv(high_label,"./outdata/network/label_c1.csv",row.names = F)
write.csv(low_label,"./outdata/network/label_c2.csv",row.names = F)

write.csv(node_high,file = "./outdata/network/network_c1.csv",row.names = F)
write.csv(node_low,file = "./outdata/network/network_c2.csv",row.names = F)


####连通度>=3 的基因做生存分析
num_high <- as.data.frame(table(data.frame(gene=c(node_high[,1],node_high[,2])))) %>% filter(Freq>=3)
num_high <- intersect(num_high[,1],Bulk_exp$TCGA$SYMBOL)
num_low <- as.data.frame(table(data.frame(gene=c(node_low[,1],node_low[,2])))) %>% filter(Freq>=3)
num_low <- intersect(num_low[,1],Bulk_exp$TCGA$SYMBOL)
###C1
net_hub_C1 <-list()
net_hub_C1$TCGA <- do.call(rbind,lapply(num_high, function(x){
  Gene_surplot(Bulk_exp$TCGA,Bulk_data$TCGA$pdata,time = "OS.time",event = "OS",gene = x,
               cutoff_method = "best",P_HR = T)
}))
net_hub_C1$ICGC <- do.call(rbind,lapply(num_high, function(x){
  Gene_surplot(Bulk_exp$ICGC,Bulk_data$ICGC$pdata,time = "OS.time",event = "OS",gene = x,
               cutoff_method = "best",P_HR = T)
}))
net_hub_C1$"In-house" <- do.call(rbind,lapply(num_high, function(x){
  Gene_surplot(Bulk_exp$"In-house",Bulk_data$"In-house"$pdata,time = "OS.time",event = "OS",gene = x,
               cutoff_method = "best",P_HR = T)
}))
##C2
net_hub_C2 <-list()
net_hub_C2$TCGA <- do.call(rbind,lapply(num_low, function(x){
  Gene_surplot(Bulk_exp$TCGA,Bulk_data$TCGA$pdata,time = "OS.time",event = "OS",gene = x,
               cutoff_method = "best",P_HR = T)
}))
net_hub_C2$ICGC <- do.call(rbind,lapply(num_low, function(x){
  Gene_surplot(Bulk_exp$ICGC,Bulk_data$ICGC$pdata,time = "OS.time",event = "OS",gene = x,
               cutoff_method = "best",P_HR = T)
}))
net_hub_C2$"In-house" <- do.call(rbind,lapply(num_low, function(x){
  Gene_surplot(Bulk_exp$"In-house",Bulk_data$"In-house"$pdata,time = "OS.time",event = "OS",gene = x,
               cutoff_method = "best",P_HR = T)
}))
#
Bulk_exp1<-Bulk_exp[c("TCGA","ICGC","In-house")]
Bulk_data1<-Bulk_data[c("TCGA","ICGC","In-house")]
hub_de_c1<-mapply(function(z,y){
  data<-egg_us(z)[num_high,]
  apply(data,1,function(x){
    label=y$predict$Subtype
    pvalue <- wilcox.test(x[which(label=="C1")],x[which(label=="C2")])$p.value
    C1median <- median(x[which(label=="C1")])
    C2median <- median(x[which(label=="C2")])
    df <- C1median-C2median
    out=df/abs(df) * pvalue
  })
},Bulk_exp1,Bulk_data1,SIMPLIFY = F)
hub_de_c1<-do.call(cbind,lapply(hub_de_c1, function(x){x}))

hub_de_c2<-mapply(function(z,y){
  data<-egg_us(z)[num_low,]
  apply(data,1,function(x){
    label=y$predict$Subtype
    pvalue <- wilcox.test(x[which(label=="C1")],x[which(label=="C2")])$p.value
    C1median <- median(x[which(label=="C1")])
    C2median <- median(x[which(label=="C2")])
    df <- C1median-C2median
    out=df/abs(df) * pvalue
  })
},Bulk_exp1,Bulk_data1,SIMPLIFY = F)
hub_de_c2<-do.call(cbind,lapply(hub_de_c2, function(x){x}))


###kegg
library(clusterProfiler)
library(org.Hs.eg.db)
cell_name<-c("CD19..B.cells","CD4..T.cells","CD8..T.cells","Eosinophils","Macrophage.Monocyte.derived","Mast.cells","NK.cells","Neutrophils","Dendritic_cells")
node_low_2<-node_low
node_low_2_unique<-unique(c(node_low_2$From,node_low_2$To))
node_low_2_unique<-setdiff(node_low_2_unique,cell_name)
node_low_2_unique<-bitr(node_low_2_unique,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")[,2]
node_low_2_unique<-as.character(node_low_2_unique)
node_low_kegg<-enrichKEGG(node_low_2_unique,use_internal_data=TRUE)
##
node_high_2<-node_high
node_high_2_unique<-unique(c(node_high_2$From,node_high_2$To))
node_high_2_unique<-setdiff(node_high_2_unique,cell_name)
node_high_2_unique<-bitr(node_high_2_unique,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")[,2]
node_high_2_unique<-as.character(node_high_2_unique)
node_high_kegg<-enrichKEGG(node_high_2_unique,use_internal_data=TRUE)

setdiff(node_low_kegg[,2],intersect(node_high_kegg[,2],node_low_kegg[,2]))

setdiff(node_high_kegg[,2],intersect(node_high_kegg[,2],node_low_kegg[,2]))

library(data.table)
kegg_info<-fread("./data/kegg_info.txt",data.table = F)
kegg_c2<-node_low_kegg[kegg_info$`big annotion`[match(node_low_kegg[,2],kegg_info$Pathway)]!="Human Diseases",]
kegg_c1<-node_high_kegg[kegg_info$`big annotion`[match(node_high_kegg[,2],kegg_info$Pathway)]!="Human Diseases",]
data.frame(name=kegg_c1[,2],GeneRatio=kegg_c1$Count/72,type=1)
kegg<-list(C1=data.frame(name=kegg_c1[,2],GeneRatio=kegg_c1$Count/72*1),
           C2=data.frame(name=kegg_c2[,2],GeneRatio=kegg_c2$Count/67*(-1))
)
kegg <- Reduce(function(x,y){full_join(x,y,by="name")},kegg)
kegg<-egg_us(kegg)
colnames(kegg)<-c("C1","C2")
kegg1 <- kegg[is.na(kegg[,1]) | is.na(kegg[,2]),]
kegg2 <- kegg[!is.na(kegg[,1]) & !is.na(kegg[,2]),]
kegg3<-rbind(kegg2,kegg1)
library(pheatmap)
p1<-pheatmap(t(kegg3),cluster_cols=F,cluster_rows=F,breaks = -c(seq(0.42,0,-0.01),seq(-0.01,-0.25,-0.01)),
         color =   c(colorRampPalette(colors = c("#0F2C67","#7989AA"))(42),
                     colorRampPalette(colors = c("#E37E7E","#CD1818"))(25)),
         na_col = "white",border_color=NA) 
pdf("./plot/kegg_network.pdf",height = 4.1,width = 5.5)
print(p1)
dev.off()

# 4.Kegg and Go ------------------------------------------------------------------

set.seed(123)
library(data.table)
library(sqldf)
library(dplyr)
library(clusterProfiler)
library(KEGG.db)
library(pheatmap)
##
library(org.Hs.eg.db)
Bulk_ID<-readRDS("./data/Bulk_ID.rds")
kegg_result <- mapply(function(x,y){
  gse_kegg(x,y$Subtype)
}, Bulk_ID,predict_data,SIMPLIFY = F)


NES <- lapply(kegg_result, function(x){
  x[,c(2,5)]
})
NES <- Reduce(function(x,y){full_join(x,y,by="Description")},NES)
colnames(NES)[-1] <- c("TCGA","ICGC","GSE14520","GSE116174","In-house")
ID<-do.call(rbind,lapply(kegg_result, function(x){
  x[,c(1,2)]
}))
ID<-unique(ID)
rownames(NES)<-ID$ID[match(NES[,1], ID$Description)]

test <- apply(NES[,-1],1,function(x){
  if(length(na.omit(x))<3){
    out=0
  }
  else{
    out=sum((as.numeric(na.omit(x)>0) -0.5)*2) / length(na.omit(x))
  }
  out
})
NES1 <- NES[test==1 | test==-1,]
NES2 <- as.matrix(NES1[,-1])
rownames(NES2) <- NES1[,1]
temp1 <- apply(NES2, 1, function(x){
  sum((as.numeric(na.omit(x)>0) -0.5)*2)
})
NES2<-NES2[order(temp1,decreasing=T),]


kegg_info<-fread("./data/kegg_info.txt",data.table = F)
main_class<-kegg_info$`big annotion`[match(rownames(NES2),kegg_info$Pathway)]
NES2 <- NES2[-which(main_class=="Human Diseases"),]
NES_type <- main_class[-which(main_class=="Human Diseases")]
NES_type <- as.data.frame(NES_type)
rownames(NES_type) <- rownames(NES2)
colnames(NES_type) <- "Mainclass"
##画图  
library(ggsci)
library(pheatmap)
ann_colors = list(
  Mainclass = c("Cellular Processes"  = pal_nejm("default")(5)[1],
                "Genetic Information Processing" = pal_nejm("default")(5)[2], 
                "Organismal Systems" = pal_nejm("default")(5)[3],
                "Environmental Information Processing" = pal_nejm("default")(5)[5],
                "Metabolism" = pal_nejm("default")(5)[4])
)
p1 <- pheatmap(NES2,cluster_cols=F,cluster_rows=F,breaks = -c(seq(3.2,0,-0.1),seq(-0.1,-2.3,-0.1)),
               color =   c(colorRampPalette(colors = c("#0F2C67","white","#CD1818"))(56)),
               na_col = "white",annotation_row = NES_type,
               annotation_colors=ann_colors) 

pdf("./plot/NES1.pdf",width=7.5,height=12)
print(p1)
dev.off()

####GO
library(dplyr)
library(GSVA)
library(clusterProfiler)
GO_immune <- fread("./data/GO_immune.txt",data.table = F)
GO_immune_set <- split(GO_immune$gene, GO_immune$term)

out <- list() 
name <- c("TCGA","ICGC","GSE14520","GSE116174","In-house")
for (i in 1:5) {
  data <- Bulk_exp[[name[i]]]
  data <- egg_us(data)
  eggtest <- gsva(data,GO_immune_set,method='ssgsea')
  temp <- data.frame(predict_data[[name[i]]]$Subtype,t(eggtest))  
  temp1 <- reshape2::melt(temp)
  temp1 <- data.frame(temp1,dataset=name[i])
  out[[name[i]]] <- temp1
}
Imm_feature <- do.call(rbind,lapply(out,function(x) x))
colnames(Imm_feature)[1]<-"subtype"
Imm_feature$dataset <- factor(Imm_feature$dataset,levels = c("TCGA","ICGC","GSE14520","GSE116174","In-house"))
Imm_feature$subtype <- factor(Imm_feature$subtype,levels = c("C1","C2"))

name <- as.character(unique(Imm_feature$variable))
name1 <- names(GO_immune_set)
plot_imm <- list() 
for (i in 1:4) {
  test<-Imm_feature[Imm_feature$variable == name[i],]
  test$variable <- as.character(test$variable)
  plot_imm[[name1[i]]] <- plot.violin(test,x=dataset,y=value,
                                      fill = subtype,
                                      palette.fill=c("#CD1818","#0F2C67"),
                                      palette.color = c("#CD1818","#0F2C67"),
                                      color = "white",
                                      fill.b = "white",
                                      color.b = subtype,
                                      group=subtype,
                                      alpha.v = 0.7,
                                      notch.b=T,
                                      width.b=0.5,
                                      add = c("boxplot","violinplot","shadow"),#,"shadow"
                                      shadow.col = c("grey40","white"),
                                      x_label="",
                                      y_label = "Score",
                                      theme = theme_bw(base_size = 16))+
    theme(axis.text.x = element_blank(),   ## 删去所有刻度标签
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 1,unit = 'cm'))+
    theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5))+
    coord_cartesian(clip = 'off')
}
plot_imm$GOBP_ACTIVATION_OF_INNATE_IMMUNE_RESPONSE
plot_imm$GOBP_ADAPTIVE_IMMUNE_RESPONSE
plot_imm$GOBP_REGULATION_OF_HUMORAL_IMMUNE_RESPONSE
plot_imm$GOBP_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS

TCGA<-gsea_egg(Bulk_exp$TCGA,predict_data$TCGA$Subtype,GO_immune)
ICGC<-gsea_egg(Bulk_exp$ICGC,predict_data$ICGC$Subtype,GO_immune)
GSE14520<-gsea_egg(Bulk_exp$GSE14520,predict_data$GSE14520$Subtype,GO_immune,log=T)
GSE116174<-gsea_egg(Bulk_exp$GSE116174,predict_data$GSE116174$Subtype,GO_immune,log=T)
Inhouse<-gsea_egg(Bulk_exp$"In-house",predict_data$"In-house"$Subtype,GO_immune)
Go_result<-list(TCGA=TCGA,ICGC=ICGC,GSE14520=GSE14520,GSE116174=GSE116174,
                "In-house"=Inhouse)

dat1<-TCGA@result[c("ID","NES","p.adjust")]
dat2<-ICGC@result[c("ID","NES","p.adjust")]
dat3<-GSE14520@result[c("ID","NES","p.adjust")]
dat4<-GSE116174@result[c("ID","NES","p.adjust")]
dat5<-Inhouse@result[c("ID","NES","p.adjust")]

NES <- lapply(Go_result, function(x){
  out<-x@result[c("ID","NES","p.adjust")]
  out<-out[out[,3]<0.05,-3]
})
NES <- Reduce(function(x,y){full_join(x,y,by="ID")},NES)
colnames(NES)[-1] <- c("TCGA","ICGC","GSE14520","GSE116174","In-house")
NES2<-egg_us(NES)
##plot
rownames(NES2)
library(ggsci)
library(pheatmap)
NES_type<-as.data.frame(rownames(NES2))
colnames(NES_type)<-"GO"
rownames(NES_type)<-NES_type[,1]
ann_colors = list(
  GO = c("GOBP_ACTIVATION_OF_INNATE_IMMUNE_RESPONSE" = pal_nejm("default")(5)[2],
         "GOBP_REGULATION_OF_HUMORAL_IMMUNE_RESPONSE"  = pal_nejm("default")(5)[1],
                "GOBP_ADAPTIVE_IMMUNE_RESPONSE" = pal_nejm("default")(5)[3],
                "GOBP_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS" = pal_nejm("default")(5)[4]))
  
p3<-pheatmap(NES2[c(2,1,3,4),],cluster_cols=F,cluster_rows=F,breaks = -c(seq(3.2,0,-0.1),seq(-0.1,-2.3,-0.1)),
               color =   c(colorRampPalette(colors = c("#0F2C67","white","#CD1818"))(56)),
               na_col = "white",annotation_row = NES_type,
               annotation_colors=ann_colors,show_rownames=F,border_color=NA) 
pdf("./plot/Go.pdf",height =3 ,width =8 )
print(p3)
dev.off()

# 5. therapy --------------------------------------------------------------
rm(list = ls())
library(data.table)
library(sqldf)
library(dplyr)
library(survival)
library(foreach)
library(utils)
library(survminer)
library(survMisc)
library(tibble)
library(vegan)
library(patchwork)
source("Function.R")

## 5.1 TIDE, TIS and IIS ---------------------------------------------------
##
Bulk_exp <- readRDS("E:/project/HCC_Immune_Subtype/data/Bulk_exp.rds")
Bulk_data <- readRDS("E:/project/HCC_Immune_Subtype/outdata/Bulk_data.rds")
#
########TIDE
TIDE<-list()
TIDE$TCGA <- fread("./data/TCGA_TIDE.csv",data.table = F)[,c(1,4)]
TIDE$TCGA <- TIDE$TCGA[match(Bulk_data$TCGA$pdata$sample,substring(TIDE$TCGA[,1],1,15)),]
TIDE$TCGA <- egg_us(TIDE$TCGA)
TIDE$TCGA <- data.frame(subtype=Bulk_data$TCGA$predict$Subtype,TIDE=TIDE$TCGA,dataset="TCGA")

TIDE$ICGC <- fread("./data/ICGC_TIDE.csv",data.table = F)[,c(1,4)]
TIDE$ICGC <- TIDE$ICGC[match(Bulk_data$ICGC$pdata$sample,TIDE$ICGC[,1]),]
TIDE$ICGC <- egg_us(TIDE$ICGC)
TIDE$ICGC <- data.frame(subtype=Bulk_data$ICGC$predict$Subtype,TIDE=TIDE$ICGC,dataset="ICGC")

TIDE$GSE14520 <- fread("./data/GSE14520_TIDE.csv",data.table = F)[,c(1,4)]
TIDE$GSE14520 <- TIDE$GSE14520[match(Bulk_data$GSE14520$pdata[,1],TIDE$GSE14520[,1]),]
TIDE$GSE14520 <- egg_us(TIDE$GSE14520)
TIDE$GSE14520 <- data.frame(subtype=Bulk_data$GSE14520$predict$Subtype,TIDE=TIDE$GSE14520,dataset="GSE14520")

TIDE$GSE116174 <- fread("./data/GSE116174_TIDE.csv",data.table = F)[,c(1,4)]
TIDE$GSE116174 <- TIDE$GSE116174[match(colnames(Bulk_exp$GSE116174)[-1],TIDE$GSE116174[,1]),]
TIDE$GSE116174 <- egg_us(TIDE$GSE116174)
TIDE$GSE116174 <- data.frame(subtype=Bulk_data$GSE116174$predict$Subtype,TIDE=TIDE$GSE116174,dataset="GSE116174")

TIDE$"In-house" <- fread("./data/inhouse_TIDE.csv",data.table = F)[,c(1,4)]
TIDE$"In-house" <- TIDE$"In-house"[match(colnames(Bulk_exp$"In-house")[-1],TIDE$"In-house"[,1]),]
TIDE$"In-house" <- egg_us(TIDE$"In-house")
TIDE$"In-house" <- data.frame(subtype=Bulk_data$"In-house"$predict$Subtype,TIDE=TIDE$"In-house",dataset="In-house")

TIDE_all<-do.call(rbind,lapply(TIDE,function(x) x))
TIDE_all$subtype<-factor(TIDE_all$subtype,levels = c("C1","C2"))
TIDE_all$dataset<-factor(TIDE_all$dataset,levels = c("TCGA","ICGC","GSE14520","GSE116174","In-house"))

save(TIDE_all,file = "./outdata/therapy/TIDE_all.Rdata")
#######TIS/IIS
library(GSVA)
library(data.table)
library(tidyr)
library(GSEABase)
library(sqldf)
APM_24<-freadlist("./data/APM_24IMMUNE.txt")
##计算ssGSEA分数
score <- lapply(Bulk_exp, function(x){
  x<-egg_us(x)
  out<-gsva(x,APM_24,method='ssgsea')
})

score<-lapply(score, function(x){
  x<-x %>% pheatmap:::scale_rows()
})
name=names(score)
score3<-mapply(function(x,y,z){
  APM_24IMMUNE(x,y$predict$Subtype,name = z)
  
},score,Bulk_data,name,SIMPLIFY = F)
  

TIS_IIS<-do.call(rbind,lapply(score3,function(x) x))
TIS_IIS<-TIS_IIS[TIS_IIS$variable!="APM",]
#
save(TIS_IIS,file = "./outdata/therapy/TIS_IIS.Rdata")

##
TIDE <- data.frame(subtype = TIDE_all$subtype,variable = "TIDE",value = TIDE_all$TIDE,
                   Dataset = TIDE_all$dataset)
colnames(TIS_IIS)[1]<-"subtype"
all_icp <- rbind(TIDE,TIS_IIS)
all_icp$Dataset<-factor(all_icp$Dataset,levels = c("TCGA","ICGC","GSE14520","GSE116174","In-house"))
all_icp$subtype<-factor(all_icp$subtype,levels = c("C1","C2"))
save(all_icp,file = "./outdata/therapy/all_icp.Rdata")

name <- as.character(unique(all_icp$variable))
name1 <- name
plot_imm <- list() 
for (i in 1:3) {
  test<-all_icp[all_icp$variable == name[i],]
  test$variable <- as.character(test$variable)
  plot_imm[[name1[i]]] <- plot.violin(test,x=Dataset,y=value,
                                      fill = subtype,
                                      palette.fill=c("#CD1818","#0F2C67"),
                                      palette.color = c("#CD1818","#0F2C67"),
                                      color = "white",
                                      fill.b = "white",
                                      color.b = subtype,
                                      color.p = subtype,
                                      fill.p = subtype,
                                      size.p=1,
                                      group=subtype,
                                      alpha.v = 0.3,
                                      alpha.b = 0.9,
                                      alpha.p = 0.6,
                                      notch.b=T,
                                      width.b=0.8,
                                      add = c("boxplot","point"),#,"shadow"
                                      shadow.col = c("grey40","white"),
                                      x_label="",
                                      y_label = "Score",
                                      theme = theme_bw(base_size = 16))+
    theme(#axis.text.x = element_blank(),   ## 删去所有刻度标签
          #axis.ticks.x = element_blank(),
          plot.margin = margin(t = 1,unit = 'cm'))+
    theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5))+
    coord_cartesian(clip = 'off')+
    facet_wrap(~variable,scales = "free_x")+
    theme(strip.text.x = element_text(size = 15, colour = "black")) + # 设置分面的字字体大小、颜色、背景、边框，
    theme(strip.background.x = element_rect(colour = NA))+
    geom_vline(xintercept = c(1.5,2.5,3.5,4.5), lwd = 1, colour = "black")
}
library(patchwork)
pdf("./plot/TIDE_TIS_IIS.pdf",height = 10,width=7)
plot_imm$TIDE/
  plot_imm$TIS/
  plot_imm$IIS
dev.off()



## 5.2 CTRP2 ---------------------------------------------------------------
library(data.table)
library(ggpubr)
library(tidyverse)
source("Function.R")
predict_data <- readRDS("E:/project/HCC_Immune_Subtype/outdata/predict_data.rds")
CTRP2 <- list()
CTRP2$TCGA <- fread("./data/DP_TCGA_CTRP2.csv",data.table = F)
CTRP2$ICGC <- fread("./data//DP_ICGC_CTRP2.csv",data.table = F)
CTRP2$GSE14520 <- fread("./data//DP_GSE14520_CTRP2.csv",data.table = F)
CTRP2$GSE116174 <- fread("./data//DP_GSE116174_CTRP2.csv",data.table = F)
CTRP2$"In-house" <- fread("./data//DP_inhouse_CTRP2.csv",data.table = F)
CTRP2_p <- mapply(function(x,y){
  mydata <- egg_us(x)
  test <- data.frame(subtype=y$Subtype,mydata)
  test <- reshape2::melt(test)
  pvalue <- compare_means(value~subtype,test,group.by = "variable")
  aaa<-test %>% group_by(variable,subtype) %>% summarise(median = median(value))
  temp<-reshape2::dcast(aaa, variable  ~ subtype,value.var = 'median')%>%
    mutate(df=(((C2-C1)>0) -0.5)*2,pvalue=pvalue$p * df)
},CTRP2,predict_data,SIMPLIFY = FALSE)
###
p <- lapply(CTRP2_p, function(x){x[,c(1,5)]})
p <- Reduce(function(x,y){full_join(x,y,by="variable")},p)
p <- egg_us(p)

colnames(p) <- c("TCGA","ICGC","GSE14520","GSE116174","In-house")
p_c <- data.frame(p,
                  SigP_Num=apply((((abs(p)<0.05)+1-1) *(abs(p) * 1/p)),1,sum))
p_c1 <- data.frame(p,
                   C1_sig_sensitive=apply(p<0.05 & p>0,1,sum),
                   C2_sig_sensitive=apply(p>-0.05 & p<0,1,sum),
                   SigP_Num=apply((((abs(p)<0.05)+1-1) *(abs(p) * 1/p)),1,sum))
library(openxlsx)
write.csv(p_c,file = "./outdata/therapy/p_c.csv",row.names = T)
write.xlsx(CTRP2_p,file = "./outdata/therapy/CTRP2_p.xlsx")
write.csv(p_c1,"./outdata/therapy/p_c1.csv",row.names = T)
##plot
name <- list("TCGA","ICGC","GSE14520","GSE116174","In-house")
drug_name <- c("navitoclax.sorafenib..1.1.mol.mol.","sorafenib")
test <- mapply(function(x,y,z){
  mydata <- egg_us(x)
  mydata <- data.frame(mydata)
  test <- data.frame(subtype=y$Subtype,mydata[,drug_name])
  test <- reshape2::melt(test)
  test <- data.frame(test,dataset=z)
},CTRP2,predict_data,name,SIMPLIFY = FALSE)

test1 <- do.call(rbind,lapply(test, function(x) x))
test1$dataset <- factor(test1$dataset,levels = c("TCGA","ICGC","GSE14520","GSE116174","In-house"))


plot_imm <- list() 
name <- c("sorafenib","navitoclax.sorafenib..1.1.mol.mol.")
name1 <- c("Sorafenib","Sorafenib+Navitoclax")

plot_imm <- list() 
for (i in 1:2) {
  test<-test1[test1$variable == name[i],]
  test$variable <- as.character(test$variable)
  plot_imm[[name1[i]]] <- plot.violin(test,x=dataset,y=value,
                                      fill = subtype,
                                      palette.fill=c("#CD1818","#0F2C67"),
                                      palette.color = c("#CD1818","#0F2C67"),
                                      color = "white",
                                      fill.b = "white",
                                      color.b = subtype,
                                      color.p = subtype,
                                      fill.p = subtype,
                                      size.p=1,
                                      group=subtype,
                                      alpha.v = 0.3,
                                      alpha.b = 0.9,
                                      alpha.p = 0.6,
                                      notch.b=T,
                                      width.b=0.8,
                                      add = c("boxplot","point"),#,"shadow"
                                      shadow.col = c("grey40","white"),
                                      x_label="",
                                      y_label = "Score",
                                      theme = theme_bw(base_size = 16))+
    theme(#axis.text.x = element_blank(),   ## 删去所有刻度标签
      #axis.ticks.x = element_blank(),
      plot.margin = margin(t = 1,unit = 'cm'))+
    theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5))+
    coord_cartesian(clip = 'off')+
    facet_wrap(~variable,scales = "free_x")+
    theme(strip.text.x = element_text(size = 15, colour = "black")) + # 设置分面的字字体大小、颜色、背景、边框，
    theme(strip.background.x = element_rect(colour = NA))+
    geom_vline(xintercept = c(1.5,2.5,3.5,4.5), lwd = 1, colour = "black")
}
pdf(file = "./plot/Sorafenib.pdf",height = 8,width = 7)
plot_imm$Sorafenib/
  plot_imm$`Sorafenib+Navitoclax`
dev.off()


## 5.3 Immunotherapy -------------------------------------------------------------
{
  ###HCC
  GSE140901 <- readRDS("./data/GSE140901_24_HCC.rds")
  GSE93647 <- readRDS("./data/GSE93647.rds")
  ###Melanoma
  GSE35640 <- readRDS("./data/GSE35640_56_M.rds")
  GSE100797 <- readRDS("./data/GSE100797_25_M.rds")
  PRJEB23709 <- readRDS("./data/PRJEB23709_73_M.rds")
  GSE91061 <- readRDS("./data/bam_51_m.rds")
  
  saveRDS(GSE140901,file = "GSE140901_24_HCC.rds")
  ###
  library(IPS29GPS)
  GSE140901$predict <-predict_29pairs(GSE140901$exp,Pairs)
  GSE93647_predict<-predict_29pairs(as.data.frame(GSE93647),Pairs)
  table(GSE93647_predict$Subtype)
  #
  GSE35640$predict <-predict_29pairs(GSE35640$exp,Pairs)
  GSE100797$predict <-predict_29pairs(GSE100797$exp,Pairs)
  PRJEB23709$predict <-predict_29pairs(PRJEB23709$exp,Pairs)
  GSE91061$predict <-predict_29pairs(GSE91061$exp,Pairs)
  ##
  GSE140901_r<-data.frame(subtype=GSE140901$predict$Subtype,GSE140901$pdata[,c("R_NR","os_time","os_event")],dataset="GSE140901")
  
  ##Me
  GSE35640_r<-data.frame(subtype=GSE35640$predict$Subtype,dataset="GSE35640",GSE35640$pdata[,c("clinical_benefit","type")],os_time=NA,os_event=NA)
  GSE100797_r<-data.frame(subtype=GSE100797$predict$Subtype,dataset="GSE100797",GSE100797$pdata[,c("clinical_benefit","type","os_time","os_event")])
  PRJEB23709_r<-data.frame(subtype=PRJEB23709$predict$Subtype,dataset="PRJEB23709",PRJEB23709$pdata[,c("clinical_benefit","type","os_time","os_event")])
  GSE91061_r<-data.frame(subtype=GSE91061$predict$Subtype,dataset="GSE91061",GSE91061$pdata[,c("clinical_benefit","type","os_time","os_event")])
  #
  Me <- rbind(GSE35640_r,GSE100797_r,PRJEB23709_r,GSE91061_r)
  Me$type[Me$type=="aCTLA4(after)+aPD1"]<-"aCTLA4+aPD1"
  
  ###
  library(data.table)
  library(survival)
  library(survminer)
  library(survMisc)
  library(patchwork)
  
  ###GSE140901 result
  {
    HCCbar <- as.data.frame(table(GSE140901_r[,c(1,2,5)]))
    HCCbar <- ddply(HCCbar,.(dataset,subtype),transform,percent=Freq/sum(Freq)*100) 
    HCCbar$label = paste0(sprintf("%.1f", HCCbar$percent), "%")
    ggplot(data = HCCbar,aes(x = subtype,y = percent,fill=R_NR)) +
      geom_bar(stat = 'identity',
               position = position_stack(),
      ) +
      scale_y_continuous(labels = scales::percent_format(scale = 1))+ 
      geom_text(aes(label=label),position = position_stack(vjust = 0.5),color="white",
                size=5)+
      theme_bw() +
      facet_wrap(~dataset,ncol = 5)+
      # 自己定义颜色
      scale_fill_manual(values =
                          c('R' = '#0072B5FF','NR' = '#BC3C29FF')) +
      theme_bw(base_size = 14,
               base_line_size = 1,
               base_rect_size = 2) +
      theme(panel.grid = element_blank(),
            axis.text = element_text(size = 14),
            axis.text.x  = element_text(colour = c("#CD1818","#0F2C67"),size = 14),
            legend.title = element_blank())
    imm_surplot(GSE140901_r,dataset = "GSE140901",type=NULL,title = "GSE140901",class="dataset")$plot
  }
  
  
  ##ME result
  {
    library(ggplot2)
    library(plyr)
    library(ggh4x)
    library(ggsci)
    library(patchwork)
    dataset <- as.data.frame(table(Me[,1:3]))
    dataset <- dataset[dataset$Freq!=0,]
    dataset <- ddply(dataset,.(dataset,subtype),transform,percent=Freq/sum(Freq)*100) 
    dataset$label = paste0(sprintf("%.1f", dataset$percent), "%")
    dataset$dataset<-factor(dataset$dataset,levels =c("PRJEB23709",
                                                      "GSE91061","GSE100797","GSE35640"))
    Type <- as.data.frame(table(Me[,c(1,3,4)]))
    Type <- Type[Type$Freq!=0,]
    Type <- ddply(Type,.(type,subtype),transform,percent=Freq/sum(Freq)*100) 
    Type$label = paste0(sprintf("%.1f", Type$percent), "%")
    Type$type<-factor(Type$type,levels =c("aPD1","aCTLA4+aPD1","ACT","MAGE-A3"))
    
    p1<-ggplot(data = dataset,aes(x = subtype,y = percent,fill=clinical_benefit)) +
      geom_bar(stat = 'identity',
               position = position_stack(),
               #width = 0.5
      ) +
      scale_y_continuous(labels = scales::percent_format(scale = 1))+ 
      geom_text(aes(label=label),position = position_stack(vjust = 0.5),color="white")+
      theme_bw() +
      facet_wrap(~dataset,ncol = 6)+
      #facet_nested(~dataset,scales="free",axes = "all")+
      #facet_wrap2(vars(dataset),ncol = 3,axes = "all")+
      scale_fill_manual(values =
                          c('Yes' = '#0072B5FF','No' = '#BC3C29FF')) +
      theme_bw(base_size = 14,
               base_line_size = 1,
               base_rect_size = 2) +
      theme(panel.grid = element_blank(),
            axis.text = element_text(size = 14),
            legend.title = element_blank(),
            axis.text.x  = element_text(colour = c("#CD1818","#0F2C67"),size = 14))
    p2<-ggplot(data = Type,aes(x = subtype,y = percent,fill=clinical_benefit)) +
      geom_bar(stat = 'identity',
               position = position_stack(),
               # width = 0.5
      ) +
      scale_y_continuous(labels = scales::percent_format(scale = 1))+ 
      geom_text(aes(label=label),position = position_stack(vjust = 0.5),color="white")+
      theme_bw() +
      facet_wrap(~type,ncol = 6)+
      scale_fill_manual(values =
                          c('Yes' = '#0072B5FF','No' = '#BC3C29FF')) +
      theme_bw(base_size = 14,
               base_line_size = 1,
               base_rect_size = 2) +
      theme(panel.grid = element_blank(),
            axis.text = element_text(size = 14),
            legend.title = element_blank(),
            axis.text.x  = element_text(colour = c("#CD1818","#0F2C67"),size = 14))
    
    imm_surplot(Me,dataset = "PRJEB23709",type=NULL,title = "PRJEB23709",class="dataset")$plot+
      imm_surplot(Me,dataset = NULL,type="aPD1",title = "aPD1",class="type")$plot+
      imm_surplot(Me,dataset = "GSE91061",type=NULL,title = "GSE91061",class="dataset")$plot+
      imm_surplot(Me,dataset = NULL,type="aCTLA4+aPD1",title = "aCTLA4+aPD1",class="type")$plot+
      imm_surplot(Me,dataset = "GSE100797",type=NULL,title = "GSE100797",class="dataset")$plot+
      imm_surplot(Me,dataset = NULL,type="ACT",title = "ACT",class="type")$plot+
      plot_layout(ncol = 2)
    
  }
}

## 5.4 Sorafenib GSE109211----------------------------------------------------------------------
{
  GSE109211 <- readRDS("./data/GSE109211.rds")
  GSE109211$predict <- predict_29pairs(GSE109211$exp,Pairs)
  out <- data.frame(GSE109211$predict,GSE109211$pdata[match(GSE109211$predict$Sample,GSE109211$pdata$`!Sample_description`),])
  out <- out[out$type=="Sor",]
  out1 <- data.frame(out[,c("Subtype","response")],dataset="GSE109211")
  df1<-as.data.frame(table(out1))
  library(plyr)
  df1<- ddply(df1,.(dataset,Subtype),transform,percent=Freq/sum(Freq)*100)
  library(ggsci)
  library(ggh4x)
  ggplot(data=out1, mapping=aes(x="response",fill=response))+
    geom_bar(stat="count",width=1,position = position_fill())+
    geom_text(stat="count",label = scales::percent(df1$percent/100), 
              size=8, position=position_fill(vjust = 0.5),colour="white")+
    coord_polar("y", start=0)+
    theme_bw(base_size = 15)+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=15, face="bold"),
      legend.position = "top"
    )+
    theme(strip.text.x = element_text(size = 15),
          strip.placement = "outside")+
    theme(strip.background.x = element_rect(fill = "grey85", colour = "grey85"),
          panel.border=element_rect(colour ="grey85"))+
    scale_fill_nejm()+
    labs(color = "Sorafenib",fill= "Sorafenib")+
    facet_nested(~ dataset+Subtype)
}



# 6. Mutation & CNV ----------------------------------------------------------------
## mutation
{
  library(ggplotify)
  library(ggstatsplot)
  library(maftools)
  C1.maf = read.maf(maf = "./data/C1.maf")
  C2.maf = read.maf(maf = "./data/C2.maf")
  ####forestplot
  pt.vs.rt <- mafCompare(m1 = C2.maf, m2 = C1.maf, m1Name = 'C2', m2Name = 'C1'
                         , minMut = 5,useCNV=FALSE)
  pt.vs.rt[["results"]]<-pt.vs.rt[["results"]][((pt.vs.rt[["results"]]$C2>9) +(pt.vs.rt[["results"]]$C1>6))>0,]
  forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05)
  ###waterfall plot
  oncoplot(maf = C1.maf, top = 10)
  oncoplot(maf = C2.maf, top = 10)
  ###coOncoplot
  gene=c(pt.vs.rt$results$Hugo_Symbol[1:9])
  coOncoplot(m1 = C1.maf, genes = gene,m2 = C2.maf, m1Name = 'C1', m2Name = 'C2',
             removeNonMutated = TRUE)
  
  
  #####co-mutation
  temp<-C1.maf@data[,c("Hugo_Symbol","Tumor_Sample_Barcode")]
  temp1<-dcast(temp,Hugo_Symbol~Tumor_Sample_Barcode)
  temp1<-as.data.frame(temp1)
  genehigh<-c(pt.vs.rt$results$Hugo_Symbol[c(1:9)])
  temp2<-temp1[match(genehigh,temp1[,1]),]
  temp2<-na.omit(temp2)
  temp3<-as.data.frame(apply(temp2[,-1],2,sum))
  co_high<-data.frame(temp3,subtype="C1")
  colnames(co_high)<-c("status","subtype")
  co_high[which(co_high$status>2),1]<-2
  
  temp<-C2.maf@data[,c("Hugo_Symbol","Tumor_Sample_Barcode")]
  temp1<-dcast(temp,Hugo_Symbol~Tumor_Sample_Barcode)
  temp1<-as.data.frame(temp1)
  genehigh<-c(pt.vs.rt$results$Hugo_Symbol[c(1:9)])
  temp2<-temp1[match(genehigh,temp1[,1]),]
  temp2<-na.omit(temp2)
  temp3<-as.data.frame(apply(temp2[,-1],2,sum))
  co_low<-data.frame(temp3,subtype="C2")
  colnames(co_low)<-c("status","subtype")
  co_low[which(co_low$status>2),1]<-2
  
  co_all<-rbind(co_low,co_high)
  ggbarstats(co_all,x=status,y=subtype)
  
}

## CNV burden
{
  load("./data/cnv_burden.Rdata")
  library(ggh4x)
  mydata <- reshape2::melt(cnv_burden[,-3])
  
  p1 <- ggplot(mydata,aes(x=subtype,value,fill=subtype))+
    geom_violin(color=NA,position = position_dodge(1))+
    geom_boxplot(width=0.25,outlier.shape = NA,position = position_dodge(1),
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
    facet_nested(variable~.,scales = "free")+
    theme(strip.placement = "outside")+
    theme(strip.background.y = element_rect(fill = "grey85", colour = "grey85"),
          panel.border=element_rect(colour ="grey70"),
          strip.text.y = element_text(size = 14),
          legend.position = "top")
  #show plot
  p1
}


# 7. SPP1-bulk ------------------------------------------------------------
library(data.table)
library(sqldf)
library(dplyr)
library(survival)
library(foreach)
library(utils)
library(survminer)
library(survMisc)
library(tibble)
library(vegan)
library(patchwork)
library(openxlsx)
source("Function.R")
##SPP1
{
  # SPP1-related genes survival analysis
  Bulk_exp <- readRDS("E:/project/HCC_Immune_Subtype/data/Bulk_exp.rds")
  Bulk_data <- readRDS("E:/project/HCC_Immune_Subtype/outdata/Bulk_data.rds")
  TCGA <- Bulk_exp$TCGA
  p1 <- Gene_surplot(TCGA,pdata=Bulk_data$TCGA$pdata,time = "OS.time",event = "OS",gene = "SPP1",
                     cutoff_method="best",
                     legend.coord=c(0.2,0.3),count = F)
  p2 <- Gene_surplot(TCGA,pdata=Bulk_data$TCGA$pdata,time = "OS.time",event = "OS",gene = "CD44",
                     cutoff_method="best",
                     legend.coord=c(0.2,0.3),count = F)
  p3 <- Gene_surplot(TCGA,pdata=Bulk_data$TCGA$pdata,time = "OS.time",event = "OS",gene = "ITGB1",
                     cutoff_method="best",
                     legend.coord=c(0.2,0.3),count = F)
  p4 <- Gene_surplot(TCGA,pdata=Bulk_data$TCGA$pdata,time = "OS.time",event = "OS",gene = "ITGA5",
                     cutoff_method="best",
                     legend.coord=c(0.2,0.3),count = F)
  p1$plot+p2$plot+p3$plot+p4$plot
  
  Inhouse <- Bulk_exp$"In-house"
  p5 <- Gene_surplot(Inhouse,pdata=Bulk_data$"In-house"$pdata,time = "OS.time",event = "OS",gene = "SPP1",
                     cutoff_method="best",
                     legend.coord=c(0.2,0.3),count = F)
  p6 <- Gene_surplot(Inhouse,pdata=Bulk_data$"In-house"$pdata,time = "OS.time",event = "OS",gene = "CD44",
                     cutoff_method="best",
                     legend.coord=c(0.2,0.3),count = F)
  p7 <- Gene_surplot(Inhouse,pdata=Bulk_data$"In-house"$pdata,time = "OS.time",event = "OS",gene = "ITGB1",
                     cutoff_method="best",
                     legend.coord=c(0.2,0.3),count = F)
  p8 <- Gene_surplot(Inhouse,pdata=Bulk_data$"In-house"$pdata,time = "OS.time",event = "OS",gene = "ITGA5",
                     cutoff_method="best",
                     legend.coord=c(0.2,0.3),count = F)
  p5$plot+p6$plot+p7$plot+p8$plot
  pdf("./plot/SPP1s_sur.pdf",height = 14,width = 7)
  p1$plot/p2$plot/p3$plot/p4$plot/p5$plot/p6$plot/p7$plot/p8$plot+plot_layout(ncol = 2)
  dev.off()
  
  ##expression of SPP1-related genes
  library(tidyverse)
  library(ggpubr)
  library(ggsci)
  library(introdataviz)
  exp<-log2(t(egg_us(TCGA[match(c("SPP1","CD44","ITGB1","ITGA5"),TCGA[,1]),]))+1)
  df <- data.frame(Subtype=Bulk_data$TCGA$predict$Subtype,
                   exp)
  colnames(df)[1] <- "Subtype"
  df <- reshape2::melt(df)
  # plot
  p1<-ggplot(df,aes(x = variable,y = value,fill = Subtype)) +
    # split violin
    geom_split_violin(alpha = 1, trim = F,color = "NA",width = 1) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2),color ="white") +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2),color ="white") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90,color = 'black',hjust = 1),
          legend.position = 'top') +
    scale_fill_manual(values=c("C1"="#CD1818","C2"="#0F2C67"))+
    scale_colour_manual(values=c("C1"="#CD1818","C2"="#0F2C67"))+
    # 添加显著性标记
    stat_compare_means(aes(group=Subtype),
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1,1),
                                        symbols = c("***", "**", "*","+", "NS")),label = "p.signif",
                       size = 6)
  p1
}

##Half violin
##expression of SPP1-related genes
library(tidyverse)
library(ggpubr)
library(ggsci)
library(introdataviz)
exp<-log2(t(egg_us(TCGA[match(c("SPP1","CD44","ITGB1","ITGA5"),TCGA[,1]),]))+1)
df <- data.frame(Subtype=Bulk_data$TCGA$predict$Subtype,
                 exp)
colnames(df)[1] <- "Subtype"
df <- reshape2::melt(df)
# plot
p1<-ggplot(df,aes(x = variable,y = value,fill = Subtype)) +
  # split violin
  geom_split_violin(alpha = 1, trim = F,color = "NA",width = 1) +
  # mean point
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2),color ="white") +
  # errorbar
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2),color ="white") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90,color = 'black',hjust = 1),
        legend.position = 'top') +
  scale_fill_manual(values=c("C1"="#CD1818","C2"="#0F2C67"))+
  scale_colour_manual(values=c("C1"="#CD1818","C2"="#0F2C67"))+
  # 添加显著性标记
  stat_compare_means(aes(group=Subtype),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1,1),
                                      symbols = c("***", "**", "*","+", "NS")),label = "p.signif",
                     size = 6)

exp<-log2(t(egg_us(Inhouse[match(c("SPP1","CD44","ITGB1","ITGA5"),Inhouse[,1]),]))+1)
df <- data.frame(Subtype=Bulk_data$"In-house"$predict$Subtype,
                 exp)
colnames(df)[1] <- "Subtype"
df <- reshape2::melt(df)
# plot
p2<-ggplot(df,aes(x = variable,y = value,fill = Subtype)) +
  # split violin
  geom_split_violin(alpha = 1, trim = F,color = "NA",width = 1) +
  # mean point
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2),color ="white") +
  # errorbar
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2),color ="white") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90,color = 'black',hjust = 1),
        legend.position = 'top') +
  scale_fill_manual(values=c("C1"="#CD1818","C2"="#0F2C67"))+
  scale_colour_manual(values=c("C1"="#CD1818","C2"="#0F2C67"))+
  # 添加显著性标记
  stat_compare_means(aes(group=Subtype),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1,1),
                                      symbols = c("***", "**", "*","+", "NS")),label = "p.signif",
                     size = 6)
pdf("./plot/SPP1s_de.pdf",height =7 ,width = 4)
p1/p2
dev.off()



###
library(GSVA)
SPP1_geneset <- freadlist("./data/SPP1_geneset.txt")
Geneset_score<-lapply(Bulk_exp, function(x){
  data <- egg_us(x)
  out<-gsva(data,SPP1_geneset,method='ssgsea')
})
##
spp1_TCGA <- egg_us(Bulk_exp$TCGA)[c("SPP1","CD44","CD68","ITGB1","ITGA5"),]
spp1_inhouse <- egg_us(Bulk_exp$`In-house`)[c("SPP1","CD44","CD68","ITGB1","ITGA5"),]
spp1_TCGA <- as.data.frame(cbind(t(spp1_TCGA),t(Geneset_score$TCGA)))
spp1_inhouse <- as.data.frame(cbind(t(spp1_inhouse),t(Geneset_score$`In-house`)))

rho1 <- apply(spp1_TCGA[,c(1,2,4,5,3,11)], 2, function(x){
  a<-cor.test(spp1_TCGA$SPP1,x,method = "spearman")[["estimate"]][["rho"]]
  b<-cor.test(spp1_TCGA$CD44,x,method = "spearman")[["estimate"]][["rho"]]
  cbind(a,b)
})
rho2 <- apply(spp1_inhouse[,c(1,2,4,5,3,11)], 2, function(x){
  a<-cor.test(spp1_inhouse$SPP1,x,method = "spearman")[["estimate"]][["rho"]]
  b<-cor.test(spp1_inhouse$CD44,x,method = "spearman")[["estimate"]][["rho"]]
  cbind(a,b)
})
temp<-rbind(rho1,rho2)
rownames(temp)<-c("TCGA-SPP1","TCGA-CD44","In-hosue-SPP1","In-hosue-CD44")
p1<-pheatmap::pheatmap(temp,cluster_cols=F,cluster_rows=F,
                       color = c(colorRampPalette(colors = c("#F5D0D0","#CD1818"))(20)),
                       border_color="white",gaps_row = 2,
                       display_numbers=T,number_color="white",fontsize_number=13)

pdf("./plot/SPP1_heatmap.pdf",height = 4.5,width = 5)
print(p1)
dev.off()

rho3 <- apply(spp1_TCGA[,c(6:10)], 2, function(x){
  a<-cor.test(spp1_TCGA$SPP1,x,method = "spearman")[["estimate"]][["rho"]]
  b<-cor.test(spp1_TCGA$CD44,x,method = "spearman")[["estimate"]][["rho"]]
  d<-cor.test(spp1_TCGA$ITGB1,x,method = "spearman")[["estimate"]][["rho"]]
  e<-cor.test(spp1_TCGA$ITGA5,x,method = "spearman")[["estimate"]][["rho"]]
  cbind(a,b,d,e)
})
rho4 <- apply(spp1_inhouse[,c(6:10)], 2, function(x){
  a<-cor.test(spp1_inhouse$SPP1,x,method = "spearman")[["estimate"]][["rho"]]
  b<-cor.test(spp1_inhouse$CD44,x,method = "spearman")[["estimate"]][["rho"]]
  d<-cor.test(spp1_inhouse$ITGB1,x,method = "spearman")[["estimate"]][["rho"]]
  e<-cor.test(spp1_inhouse$ITGA5,x,method = "spearman")[["estimate"]][["rho"]]
  cbind(a,b,d,e)
})

rho5 <- apply(spp1_TCGA[,c(6:10)], 2, function(x){
  a<-cor.test(spp1_TCGA$SPP1,x,method = "spearman")$p.value
  b<-cor.test(spp1_TCGA$CD44,x,method = "spearman")$p.value
  d<-cor.test(spp1_TCGA$ITGB1,x,method = "spearman")$p.value
  e<-cor.test(spp1_TCGA$ITGA5,x,method = "spearman")$p.value
  cbind(a,b,d,e)
})
rho6 <- apply(spp1_inhouse[,c(6:10)], 2, function(x){
  a<-cor.test(spp1_inhouse$SPP1,x,method = "spearman")$p.value
  b<-cor.test(spp1_inhouse$CD44,x,method = "spearman")$p.value
  d<-cor.test(spp1_inhouse$ITGB1,x,method = "spearman")$p.value
  e<-cor.test(spp1_inhouse$ITGA5,x,method = "spearman")$p.value
  cbind(a,b,d,e)
})

rho5[,]<-apply(rho5,2,function(x){
  x[x>=0.05]<-NA
  x
} )
rho3<-rho3*((rho5+1)/(rho5+1))

rho6[,]<-apply(rho6,2,function(x){
  x[x>=0.05]<-NA
  x
} )
rho4<-rho4*((rho6+1)/(rho6+1))

out<-rbind(rho3,rho4)
rownames(out)<-c("SPP1 (TCGA)","CD44 (TCGA)","ITGB1 (TCGA)","ITGA5 (TCGA)",
                 "SPP1 (In-house)","CD44 (In-house)","ITGB1 (In-house)","ITGA5 (In-house)")
p2<-pheatmap::pheatmap(out,cluster_cols=F,cluster_rows=F,
                       color = c(colorRampPalette(colors = c("#F5D0D0","#CD1818"))(20)),
                       border_color="white",gaps_row = 4,
                       display_numbers=T,number_color="white",fontsize_number=13,na_col = "white")
pdf("./plot/SPP1_heatmap_feature.pdf",height = 6,width = 5)
print(p2)
dev.off()




