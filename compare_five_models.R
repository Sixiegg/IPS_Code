library(data.table)
library(survminer)
library(survival)


# 1. method compare (the risk thresholds from the respective dataset) -------------------------------------------------------
Bulk_exp <- readRDS("./data/Bulk_exp.rds")
Bulk_data <- readRDS("E:/project/HCC_Immune_Subtype/outdata/Bulk_data.rds")
predict_data <- readRDS("E:/project/HCC_Immune_Subtype/outdata/predict_data.rds")
data <- Bulk_exp$GSE116174
rownames(data) <- data[,1]
data <- data[,-1]
data <- log(data+1)

cal_method_score <- function(data){
  if(max(data)>20){
    data <- log(data+1)
  }
  # 0.011 * RACGAP1 - 0.024 * HAO2 - 0.055* OGDHL + 0.122*ZWINT - 0.069*CFHR3 - 0.044*CYP2C9+0.07*SFN-0.005*SPP2
  Liu_2022 <- data.frame(gene=c("RACGAP1", "HAO2", "OGDHL", "ZWINT", "CFHR3", "CYP2C9", "SFN", "SPP2"),
                         value=c(0.011, - 0.024, - 0.055, 0.122, -0.069,- 0.044, 0.07, -0.005))
  inter_gene <- intersect(Liu_2022$gene, rownames(data))
  Liu_2022_score <- data[match(inter_gene, rownames(data)),] * Liu_2022$value[match(inter_gene, Liu_2022$gene)]
  Liu_2022_score <- apply(Liu_2022_score, 2, sum)
  
  #RiskScore = 0.257 * Expression CSTB + 0.263 * Expression TALDO1 + 0.313 * Expression CLTA
  Lu_2022 <- data.frame(gene=c("CSTB", "TALDO1", "CLTA"),
                        value=c(0.257, 0.263, 0.313))
  inter_gene <- intersect(Lu_2022$gene, rownames(data))
  Lu_2022_score <- data[match(inter_gene, rownames(data)),] * Lu_2022$value[match(inter_gene, Lu_2022$gene)]
  Lu_2022_score <- apply(Lu_2022_score, 2, sum)
  
  #Risk score= -0.19231CFH+0.5966RAP1A+0.22919ENO1-0.19329GP1BA+0.12790*SLC2A1
  Liu_2023 <- data.frame(gene=c("CFH", "RAP1A", "ENO1","GP1BA","SLC2A1"),
                         value=c(-0.19231, 0.5966, 0.22919, - 0.19329, 0.12790))
  inter_gene <- intersect(Liu_2023$gene, rownames(data))
  Liu_2023_score <- data[match(inter_gene, rownames(data)),] * Liu_2023$value[match(inter_gene, Liu_2023$gene)]
  Liu_2023_score <- apply(Liu_2023_score, 2, sum)
  
  #RiskScore = -0.088ADAMTSL2 - 0.121SLCO2A1 - 0.217CD4 + 0.249GCN1 + 0.345HMGXB3 + 0.271LUC7L3
  Yu_2022 <- data.frame(gene=c("ADAMTSL2", "SLCO2A1", "CD4","GCN1","HMGXB3", "LUC7L3"),
                        value=c(-0.088, - 0.121, - 0.217, 0.249, 0.345,0.271))
  inter_gene <- intersect(Yu_2022$gene, rownames(data))
  Yu_2022_score <- data[match(inter_gene, rownames(data)),] * Yu_2022$value[match(inter_gene, Yu_2022$gene)]
  Yu_2022_score <- apply(Yu_2022_score, 2, sum)
  
  
  #risk score = (0.400506 × DBF4 expression) + (0.188240 × ARG2 expression) + (0.204192 × SLC16A3 expression)
  Tang_2022 <- data.frame(gene=c("DBF4", "ARG2", "SLC16A3"),
                          value=c(0.400506, 0.188240, 0.204192))
  inter_gene <- intersect(Tang_2022$gene, rownames(data))
  Tang_2022_score <- data[match(inter_gene, rownames(data)),] * Tang_2022$value[match(inter_gene, Tang_2022$gene)]
  Tang_2022_score <- apply(Tang_2022_score, 2, sum)
  
  out <- data.frame(Liu_2022_score=Liu_2022_score, Lu_2022_score=Lu_2022_score, Liu_2023_score=Liu_2023_score,
                    Yu_2022_score=Yu_2022_score, Tang_2022_score=Tang_2022_score)
  return(out)
}
cal_score <- lapply(Bulk_exp, function(data){
  rownames(data) <- data[,1]
  data <- data[,-1]
  cal_method_score(data)
})
  
P_HR_C_other <- lapply(c("TCGA","ICGC","GSE14520", "GSE116174", "In-house"), function(x){
  temp <- cal_score[[x]]
  index_name <- colnames(temp)
  out <- lapply(index_name, function(y){
    feature_surplot(data = temp[[y]], pdata = Bulk_data[[x]]$pdata,
                    time = "OS.time", event ="OS",
                    cutoff_method = "median", title = "data_score",P_HR = T)
  })
  out <- do.call(rbind,out)
  rownames(out) <- index_name
  return(out)
})

P_HR_C_pair <- lapply(c("TCGA","ICGC","GSE14520", "GSE116174", "In-house"), function(x){
  out <- feature_surplot(data = predict_data[[x]]$Ratio, pdata = Bulk_data[[x]]$pdata,
                  time = "OS.time", event ="OS",
                  cutoff_method = "pair", title = "data_score",P_HR = T)
  rownames(out) <- 'IPS29'
  return(out)
})
P_HR_C_rbind <- mapply(function(x,y){
  out <- rbind(y,x)
}, x=P_HR_C_other, y=P_HR_C_pair,SIMPLIFY = F)

C_index <- do.call(cbind,lapply(P_HR_C_rbind,function(x) x[,"C_index"]))
colnames(C_index) <- c("TCGA","ICGC","GSE14520", "GSE116174", "In-house")
C_index <- data.frame(method=rownames(C_index), C_index)
C_index <-melt(C_index)

p_index <- do.call(cbind,lapply(P_HR_C_rbind,function(x) x[,"p.val"]))
colnames(p_index) <- c("TCGA","ICGC","GSE14520", "GSE116174", "In-house")
p_index <- data.frame(method=rownames(p_index), p_index)
p_index <-melt(p_index)

C_index$pvalue<-ifelse(p_index$value<0.05, "p < 0.05", "p > 0.05")
C_index$method <- factor(C_index$method, levels = unique(C_index$method)[6:1])


p1 <- ggplot(C_index, aes( value, method, fill =pvalue )) +
  geom_bar(stat = "identity", position =  "dodge") +
  # theme_economist(base_size = 14) + scale_fill_economist() +
  theme(axis.ticks.length =
          unit(0.5, 'cm')) +
  scale_fill_manual(values = c("p < 0.05" = "#1F78B4", "p > 0.05" = "#A6CEE3")) +
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black"))+
  labs(x="C-index")+
  facet_grid(. ~ variable)

pdf("./plot/compare/compare_cindex.pdf", width = 10,height = 2.5)
p1
dev.off()


# 2. method compare (the risk thresholds from TCGA as a benchmark) -------------------------------------------------------------
cutoff <- apply(cal_score$TCGA, 2, median)

P_HR_C_other <- lapply(c("TCGA","ICGC","GSE14520", "GSE116174", "In-house"), function(x){
  temp <- cal_score[[x]]
  index_name <- colnames(temp)
  out <- lapply(index_name, function(y){
    feature_surplot(data = temp[[y]], pdata = Bulk_data[[x]]$pdata,
                    time = "OS.time", event ="OS",
                    cutoff_method = "cutoff", median_cutoff=cutoff[y],P_HR = T)
  })
  out <- do.call(rbind,out)
  rownames(out) <- index_name
  return(out)
})

P_HR_C_pair <- lapply(c("TCGA","ICGC","GSE14520", "GSE116174", "In-house"), function(x){
  out <- feature_surplot(data = predict_data[[x]]$Ratio, pdata = Bulk_data[[x]]$pdata,
                         time = "OS.time", event ="OS",
                         cutoff_method = "pair", title = "data_score",P_HR = T)
  rownames(out) <- 'IPS29'
  return(out)
})
P_HR_C_rbind <- mapply(function(x,y){
  out <- rbind(y,x)
}, x=P_HR_C_other, y=P_HR_C_pair,SIMPLIFY = F)

C_index <- do.call(cbind,lapply(P_HR_C_rbind,function(x) x[,"C_index"]))
colnames(C_index) <- c("TCGA","ICGC","GSE14520", "GSE116174", "In-house")
C_index <- data.frame(method=rownames(C_index), C_index)
C_index <-melt(C_index)

p_index <- do.call(cbind,lapply(P_HR_C_rbind,function(x) x[,"p.val"]))
colnames(p_index) <- c("TCGA","ICGC","GSE14520", "GSE116174", "In-house")
p_index <- data.frame(method=rownames(p_index), p_index)
p_index <-melt(p_index)

C_index$pvalue<-ifelse(p_index$value<0.05, "p < 0.05", "p > 0.05")
C_index$method <- factor(C_index$method, levels = unique(C_index$method)[6:1])


p2 <- ggplot(C_index, aes( value, method, fill =pvalue )) +
  geom_bar(stat = "identity", position =  "dodge") +
  # theme_economist(base_size = 14) + scale_fill_economist() +
  theme(axis.ticks.length =
          unit(0.5, 'cm')) +
  scale_fill_manual(values = c("p < 0.05" = "#1F78B4", "p > 0.05" = "#A6CEE3")) +
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black"))+
  labs(x="C-index")+
  facet_grid(. ~ variable)
p2
pdf("./plot/compare/TCGA_cindex.pdf", width = 10,height = 2.5)
p2
dev.off()

# 3. consistency ratio (from the respective dataset vs from all samples in five datasets merged together) -------------------------------------------------------------
cal_score <- lapply(Bulk_exp, function(data){
  rownames(data) <- data[,1]
  data <- data[,-1]
  cal_method_score(data)
})

##random add 1 value
set.seed(1)
seeds <- sample(1:1000, 50)
out <- lapply(seeds, function(z){
  score <- mapply(function(data,data_initial){
    rownames(data) <- data[,1]
    data <- data[,-1]
    n <- round(ncol(data) * 0.5)
    set.seed(z)
    sampled_cols <- sample(colnames(data), n)
    data[,sampled_cols] <- data[,sampled_cols]+1
    out <- cal_method_score(data)
    set_to_binary <- function(x) {
      median_val <- median(x)
      ifelse(x > median_val, 1, 0)
    }
    label <- apply(out, 2, set_to_binary)
    initial_label <- apply(data_initial, 2, set_to_binary)
    label_check <- (label+initial_label)!=1
    ratio <- apply(label_check, 2, function(x){
      sum(x)/length(x)
    })
  },data=Bulk_exp, data_initial =cal_score, SIMPLIFY = T)
  print(z)
  score <- data.table(method= rownames(score), score)
})

out_mean_1 <- do.call(rbind, out)
out_mean_1 <- out_mean_1[,lapply(.SD,function(x){mean(x)}) , by=method]
out_mean_1 <- rbind(out_mean_1, data.frame(method="IPS29","TCGA"=1,"ICGC"=1,"GSE14520"=1, "GSE116174"=1, "In-house"=1, check.names = F))

out <- lapply(seeds, function(z){
  score <- mapply(function(data,data_initial){
    rownames(data) <- data[,1]
    data <- data[,-1]
    n <- round(ncol(data) * 0.5)
    set.seed(z)
    sampled_cols <- sample(colnames(data), n)
    data[,sampled_cols] <- data[,sampled_cols]+5
    out <- cal_method_score(data)
    set_to_binary <- function(x) {
      median_val <- median(x)
      ifelse(x > median_val, 1, 0)
    }
    label <- apply(out, 2, set_to_binary)
    initial_label <- apply(data_initial, 2, set_to_binary)
    label_check <- (label+initial_label)!=1
    ratio <- apply(label_check, 2, function(x){
      sum(x)/length(x)
    })
  },data=Bulk_exp, data_initial =cal_score, SIMPLIFY = T)
  print(z)
  score <- data.table(method= rownames(score), score)
})
out_mean_5 <- do.call(rbind, out)
out_mean_5 <- out_mean_5[,lapply(.SD,function(x){mean(x)}) , by=method]
out_mean_5 <- rbind(out_mean_5, data.frame(method="IPS29","TCGA"=1,"ICGC"=1,"GSE14520"=1, "GSE116174"=1, "In-house"=1, check.names = F))

out <- lapply(seeds, function(z){
  score <- mapply(function(data,data_initial){
    rownames(data) <- data[,1]
    data <- data[,-1]
    n <- round(ncol(data) * 0.5)
    set.seed(z)
    sampled_cols <- sample(colnames(data), n)
    data[,sampled_cols] <- data[,sampled_cols]*2
    out <- cal_method_score(data)
    set_to_binary <- function(x) {
      median_val <- median(x)
      ifelse(x > median_val, 1, 0)
    }
    label <- apply(out, 2, set_to_binary)
    initial_label <- apply(data_initial, 2, set_to_binary)
    label_check <- (label+initial_label)!=1
    ratio <- apply(label_check, 2, function(x){
      sum(x)/length(x)
    })
  },data=Bulk_exp, data_initial =cal_score, SIMPLIFY = T)
  print(z)
  score <- data.table(method= rownames(score), score)
})
out_mean_2x <- do.call(rbind, out)
out_mean_2x <- out_mean_2x[,lapply(.SD,function(x){mean(x)}) , by=method]
out_mean_2x <- rbind(out_mean_2x, data.frame(method="IPS29","TCGA"=1,"ICGC"=1,"GSE14520"=1, "GSE116174"=1, "In-house"=1, check.names = F))

out_mean_1 <- melt(out_mean_1)
out_mean_1$type <- "value +1"
out_mean_5 <- melt(out_mean_5)
out_mean_5$type <- "value +5"
out_mean_2x <- melt(out_mean_2x)
out_mean_2x$type <- "value *2"

out_all <- rbind(out_mean_1,out_mean_5,out_mean_2x)
p3 <- ggplot(out_all, aes( method, value, fill =type )) +
  geom_bar(stat = "identity", position =  "dodge") +
  # theme_economist(base_size = 14) + scale_fill_economist() +
  theme(axis.ticks.length =
          unit(0.5, 'cm')) +
  scale_fill_manual(values = c("#1F78B4", "#33A02C","#FF7F00")) +
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"))+
  labs(y="Consistency ratio")+
  facet_grid(. ~ variable)

pdf("./plot/compare/ratio.pdf", width = 12,height = 4)
p3
dev.off()


# 4. method compare (the risk thresholds from all samples in five datasets merged together) --------------------------------------------------------------
cal_score <- lapply(Bulk_exp, function(data){
  rownames(data) <- data[,1]
  data <- data[,-1]
  cal_method_score(data)
})
set_to_binary <- function(x) {
  median_val <- median(x)
  ifelse(x > median_val, 1, 0)
}

initial_label <- mapply(function(x,y ){
  label <- apply(x, 2, set_to_binary)
  label <- data.frame(label, dataset=y)
},x=cal_score, y=names(cal_score),SIMPLIFY = F)
initial_label <- do.call(rbind, initial_label)


all_bind <- do.call(rbind,cal_score)
all_bind <- apply(all_bind, 2, set_to_binary)
label <- all_bind+initial_label[,-6] !=1
label <- data.table(label, dataset=initial_label$dataset)
ratio <- label[,lapply(.SD,function(x) sum(x,na.rm=T)/length(x)),by=dataset]
ratio$IPS29 <- 1

ratio <- melt(ratio)
ratio$variable <- factor(ratio$variable, levels = unique(ratio$variable)[c(6,1:5)])
ratio$dataset <- factor(ratio$dataset, levels = unique(ratio$dataset)[1:5])
p4 <- ggplot(ratio, aes( variable, value, fill =dataset )) +
  geom_bar(stat = "identity", position =  "dodge") +
  # theme_economist(base_size = 14) + scale_fill_economist() +
  theme(axis.ticks.length =
          unit(0.5, 'cm')) +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FDBF6F")) +
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"))+
  labs(y="Consistency ratio")+
  facet_grid(. ~ dataset)

pdf("./plot/compare/all_ratio.pdf", width = 12,height = 4)
p4
dev.off()


# 4.1 cindex --------------------------------------------------------------
all_bind <- do.call(rbind,cal_score)

cutoff <- apply(all_bind, 2, median)

P_HR_C_other <- lapply(c("TCGA","ICGC","GSE14520", "GSE116174", "In-house"), function(x){
  temp <- cal_score[[x]]
  index_name <- colnames(temp)
  out <- lapply(index_name, function(y){
    feature_surplot(data = temp[[y]], pdata = Bulk_data[[x]]$pdata,
                    time = "OS.time", event ="OS",
                    cutoff_method = "cutoff", median_cutoff=cutoff[y],P_HR = T)
  })
  out <- do.call(rbind,out)
  rownames(out) <- index_name
  return(out)
})

P_HR_C_pair <- lapply(c("TCGA","ICGC","GSE14520", "GSE116174", "In-house"), function(x){
  out <- feature_surplot(data = predict_data[[x]]$Ratio, pdata = Bulk_data[[x]]$pdata,
                         time = "OS.time", event ="OS",
                         cutoff_method = "pair", title = "data_score",P_HR = T)
  rownames(out) <- 'IPS29'
  return(out)
})
P_HR_C_rbind <- mapply(function(x,y){
  out <- rbind(y,x)
}, x=P_HR_C_other, y=P_HR_C_pair,SIMPLIFY = F)

C_index <- do.call(cbind,lapply(P_HR_C_rbind,function(x) x[,"C_index"]))
colnames(C_index) <- c("TCGA","ICGC","GSE14520", "GSE116174", "In-house")
C_index <- data.frame(method=rownames(C_index), C_index)
C_index <-melt(C_index)

p_index <- do.call(cbind,lapply(P_HR_C_rbind,function(x) x[,"p.val"]))
colnames(p_index) <- c("TCGA","ICGC","GSE14520", "GSE116174", "In-house")
p_index <- data.frame(method=rownames(p_index), p_index)
p_index <-melt(p_index)

C_index$pvalue<-ifelse(p_index$value<0.05, "p < 0.05", "p > 0.05")
C_index$method <- factor(C_index$method, levels = unique(C_index$method)[6:1])


p5 <- ggplot(C_index, aes( value, method, fill =pvalue )) +
  geom_bar(stat = "identity", position =  "dodge") +
  # theme_economist(base_size = 14) + scale_fill_economist() +
  theme(axis.ticks.length =
          unit(0.5, 'cm')) +
  scale_fill_manual(values = c("p < 0.05" = "#1F78B4", "p > 0.05" = "#A6CEE3")) +
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black"))+
  labs(x="C-index")+
  facet_grid(. ~ variable)
p5
pdf("./plot/compare/all_cindex.pdf", width = 10,height = 2.5)
p5
dev.off()


# function ----------------------------------------------------------------
Surplot <- function(data_demo,type,title){
  if(type=="OS"){
    surv_TTP<-surv_fit(Surv(OS.time, OS) ~ subtype, data =data_demo)
    xlab = "Overall Survival (days)"
  }
  else if(type=="PFI"){
    surv_TTP<-surv_fit(Surv(PFI.time, PFI) ~ subtype, data =data_demo)
    xlab = "Progression Free Interval (days)"
  }
  else{
    surv_TTP<-surv_fit(Surv(RFS.time, RFS) ~ subtype, data =data_demo)
    xlab = "Relapse Free Survival (days)"
  }
  TitleNames <- "Subtype"  
  title <- title  
  legend_label <- c(paste0("C1 n=",table(data_demo$subtype)[1]),
                    paste0("C2 n=",table(data_demo$subtype)[2]))
  ggsurvplot(surv_TTP,
             pval = TRUE, 
             pval.size = 6.5,  
             pval.coord=c(3,5),
             conf.int = F, 
             fun = "pct", 
             palette =c("#CD1818","#0F2C67"),
             risk.table = F, 
             xlab = xlab,  
             font.tickslab = 16, 
             font.x = 18,font.y = 18,  
             font.subtitle = 18, 
             font.legend = 15,  
             legend = c(0.8,0.9),  
             legend.labs = legend_label,
             censor.size = 0,  
             size = 1.8,  
             axes.offset = T, 
             legend.title = TitleNames,
             main = "Survival curves",
             submain = title)
}
feature_surplot <- function(data, pdata, feature=NA,title=NA, time, event, cutoff_method="best",
                            score_method=NA,
                            palette =c("#CD1818","#0F2C67"),
                            xlab="Overall Survival (days)",
                            count=T,
                            type.sur="default",
                            P_HR=F,
                            ...){
  #data pre-processing
  if(is.na(score_method)){
    if(is.na(feature)){
      exp <- as.numeric(data)
      feature = "feature"
    }else{
      exp <- data[feature,]
      exp <- as.numeric(exp)
    }
  }else{
    score <- data
    score <- structure(score, class = c(class(score),score_method))
    exp <- score_type(data=score, feature=feature)
  }
  #title
  if(is.na(title)){
    title=feature
  }
  
  
  
  
  #
  if(mean(exp)==0){
    return(NULL)
  }
  ##sur data
  pdata <- pdata[,c(time,event)]
  sur <- data.frame(feature=exp,pdata)
  colnames(sur)[1] <- "feature"
  
  ##
  sur <- structure(sur, class = c(class(sur),cutoff_method))
  #return(sur)
  twosurdata <- sur_cutoff(sur,time,event,...)

  if(P_HR){   
    twosurdata[,3]<-factor(twosurdata[,3],levels = c("low","high"))
    if(length( unique(twosurdata[,3]) )==1 ){
      print(1)
      FORMULA <-as.formula(paste("Surv(",time,",",event,") ~ 1"))
      p.val=1
      res.cox <- coxph(FORMULA, data = twosurdata)
      HR<-NA
      C_index <-res.cox$concordance[6]
    }else{
      FORMULA <-as.formula(paste("Surv(",time,",",event,") ~ feature"))
      sdf <- survdiff(formula=FORMULA,data=twosurdata)
      p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
      res.cox <- coxph(FORMULA, data = twosurdata)
      HR<-(summary(res.cox))[["coefficients"]][2]
      C_index <- (summary(res.cox))[["concordance"]][1]
    }
    
    out <- cbind(p.val,HR,C_index) 
    rownames(out) <- "feature"
    return(out)
  }
  
  #return(twosurdata)
  #return(twosurdata)
  FORMULA <-as.formula(paste("Surv(",time,",",event,") ~ feature"))
  surv_TTP<-surv_fit(FORMULA, data =twosurdata)
  ###
  xlab = xlab
  if(count){
    legend_label <- c(paste0("High n=",table(twosurdata[["feature"]])[1]),
                      paste0("Low n=",table(twosurdata[["feature"]])[2]))
  }
  else{
    legend_label <- c('High',"Low")
  }
  
  surv_TTP <- structure(surv_TTP, class = c(class(surv_TTP),type.sur))
  p1 <- surplot(surv_TTP,palette,title,xlab,legend_label,...)
  return(p1)
}
###



##
score_type <- function(data,feature, ...){
  UseMethod("score_type")
}
score_type.ssGSEA <- function(data,feature,...){
  library(GSVA)
  gene_set <- list(feature=feature)
  data <- as.matrix(data)
  score <- gsva(data,gene_set,method='ssgsea')
  return(as.numeric(score))
}
score_type.mean <- function(data,feature,...){
  temp <- data[feature,]
  score <- apply(temp,2,function(x){
    x <- na.omit(x)
    out <- mean(x)
    return(out)
  })
  return(as.numeric(score))
}


###
sur_cutoff <- function(sur,...){
  UseMethod("sur_cutoff")
}
##
sur_cutoff.best <- function(sur,time,event, ...){
  res.cut1 <- surv_cutpoint(sur, time=time , event =event, variables = colnames(sur)[1])
  res.cat <- surv_categorize(res.cut1)
  return(res.cat)
}
##
sur_cutoff.median <- function(sur,  ...){
  cutoff <- median(sur[,1])
  out <- data.frame(sur[,-1],ifelse(sur[,1]>cutoff,"high","low"))
  colnames(out)[3] <- colnames(sur)[1]
  return(out)
}
##
sur_cutoff.mean <- function(sur,  ...){
  cutoff <- mean(sur[,1])
  out <- data.frame(sur[,-1],ifelse(sur[,1]>cutoff,"high","low"))
  colnames(out)[3] <- colnames(sur)[1]
  return(out)
}
sur_cutoff.pair <- function(sur,  ...){
  cutoff <- 0.5
  out <- data.frame(sur[,-1],ifelse(sur[,1]>=cutoff,"high","low"))
  colnames(out)[3] <- colnames(sur)[1]
  return(out)
}

sur_cutoff.cutoff <- function(sur,median_cutoff,  ...){
  print(median_cutoff)
  cutoff <- as.numeric(median_cutoff)
  out <- data.frame(sur[,-1],ifelse(sur[,1]>cutoff,"high","low"))
  colnames(out)[3] <- colnames(sur)[1]
  return(out)
}
####
surplot <- function(surv_TTP,...){
  UseMethod("surplot")
}
surplot.default <- function(surv_TTP,palette,title,xlab,legend_label,legend.coord=c(0.8,0.9),...){
  ggsurvplot(surv_TTP,
             pval = TRUE, 
             pval.size = 6.5,  
             pval.coord=c(3,5),
             conf.int = F, 
             fun = "pct", 
             palette =palette,
             risk.table = F, 
             xlab = xlab,  
             font.tickslab = 16, 
             font.x = 18,font.y = 18,  
             font.subtitle = 18, 
             font.legend = 15,  
             legend =legend.coord,  
             legend.labs = legend_label,
             censor.size = 0,  
             size = 1.8,  
             axes.offset = T,
             main = "Survival curves",
             submain = title)
}

