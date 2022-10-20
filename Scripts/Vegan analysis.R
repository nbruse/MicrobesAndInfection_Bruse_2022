# Analysis Microbiome Bruse et al.
# Author: NB 

library(vegan)
library(ggplot2)
library(xlsx)
library(dplyr)

set.seed(123)

setwd("")

# function for rda barplot variance
barplot.veg<-function(x){
  constrained_eig <- x$CCA$eig/x$tot.chi*100
  unconstrained_eig <- x$CA$eig/x$tot.chi*100
  expl_var <- c(constrained_eig, unconstrained_eig)
  barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
           las = 2, ylab = '% variation')
}
# function to get the pvalue of permutation
get.pval<-function(x){
  file = tempfile(fileext=".txt")
  sink(file)
  print(x)
  sink(file = NULL)
  p.value<-scan(file, what="character")[37]
  return(p.value)
}
# function plot with permutation significance
overview.plot<-function(x=1,y=2,l=7){
  names<-c()
  variances<-c()
  pvals<-c()
  
  for(i in colnames(env)){
    
    env.t <- env[,i,drop=F]
    env.t<-env.t[complete.cases(env.t), ,drop=F]
    com.t<-com[rownames(com) %in% rownames(env.t),]

    
    tbRDA<-rda(as.formula(paste0("com.t","~",i)), data = env.t)
    
    pval<-permutest(tbRDA,10000)
    p.temp<-get.pval(pval)
    pvals<-c(pvals,p.temp)
  
    variance<-unname(tbRDA$CCA$eig/tbRDA$tot.chi*100)
    
    names<-c(names, i)
    variances<-c(variances, variance)
  }
  
  pvals<-as.numeric(pvals)
  df <- data.frame(variances, names, pvals)
  df<-df[order(df$variances, decreasing = T),]
  volgorde<-df$names
  df$names <- factor(df$names,levels = rev(volgorde))
  
  highlight_df <- df %>% 
    filter(pvals<0.05)
  
  if(dim(highlight_df)[1] == 0){
    test<-ggplot(df, aes(x=names, y=variances, fill=names)) + 
      geom_dotplot(binaxis='y', stackdir='center', binwidth=0.03, dotsize = y) + 
      coord_flip() + 
      theme_gray() + theme(legend.position = "none",text = element_text(size =20)) +
      xlab("") + 
      ylab("Variance explained (%)") +
      ylim(0, l)
  } else {
    test<-ggplot(df, aes(x=names, y=variances, fill=names)) + 
      geom_dotplot(binaxis='y', stackdir='center', binwidth=0.03, dotsize = y) + 
      geom_dotplot(data=highlight_df, 
                 aes(x=names,y=variances, fill = '#FF0000'), 
                 color='red', dotsize = x)+
      coord_flip() + 
      theme_gray() + theme(legend.position = "none",text = element_text(size = 15)) +
      xlab("") + 
      ylab("Variance explained (%)")
  }
  
  return(test)
}
# function plot without permutation
overview.plot2<-function(y=2,lower=0,upper=6){
  names<-c()
  variances<-c()
  pvals<-c()
  
  for(i in colnames(env)){
    
    env.t <- env[,i,drop=F]
    env.t<-env.t[complete.cases(env.t), ,drop=F]
    com.t<-com[rownames(com) %in% rownames(env.t),]
    
    
    tbRDA<-rda(as.formula(paste0("com.t","~",i)), data = env.t)
    
    variance<-unname(tbRDA$CCA$eig/tbRDA$tot.chi*100)
    
    names<-c(names, i)
    variances<-c(variances, variance)
  }
  
  df <- data.frame(variances, names)
  df<-df[order(df$variances, decreasing = T),]
  df$names<-clean(df$names)
  volgorde<-df$names
  df$names <- factor(df$names,levels = rev(volgorde))
  
  
  test<-ggplot(df, aes(x=names, y=variances, fill=names)) + 
    geom_dotplot(binaxis='y', stackdir='center', binwidth=0.03, dotsize = y) + 
    coord_flip() + 
    theme(legend.position = "none",text = element_text(size =20)) +
    xlab("") + 
    ylab("Explained variance (%)") +
    ylim(lower, upper) + 
    theme_bw() + 
    theme(legend.position="none", text = element_text(size = 18))
  
  return(test)
}
# function to potentially clean names of overview.plot
clean<-function(x){
  x<-gsub("_D0_T0", "", x)
  x<-gsub("_D0", "", x)
  x<-gsub("Max.delta", "", x)
  x<-gsub("PEAK", "", x)
  
  x<-gsub("IL6", "IL-6", x)
  x<-gsub("IL10", "IL-10", x)
  x<-gsub("IP10", "IP-10", x)
  x<-gsub("MCP1", "MCP-1", x)
  x<-gsub("IL1RA", "IL-1RA", x)
  x<-gsub("IL1b", paste0("IL-1","\U03B2"), x)
  x<-gsub("MIP1a", paste0("MIP-1","\U03B1"), x)
  x<-gsub("_Tol", " - Tolerance", x)
  x<-gsub("veggy_true", "Diet", x)
  x<-gsub("PiC", "PIC", x)
  x<-gsub("GCSF", "G-CSF", x)
  x<-gsub("IL8", "IL-8", x)
  x<-gsub("HR", "Heartrate", x)
  x<-gsub("TEMP", "Temperature", x)
  x<-gsub("ABPm", "Blood pressure", x)
  x<-gsub("ZS", "Symptom score", x)
  
  x<-gsub("_", " - ", x)
  x<-gsub("^ -", "", x)
  x<-gsub(" - $", "", x)
  
  return(x)
}
# function to read out explained variance and p val 
get_me<-function(x){
  for(i in x){
    env.t <- env[,i,drop=F]
    env.t<-env.t[complete.cases(env.t), ,drop=F]
    com.t<-com[rownames(com) %in% rownames(env.t),]
    
    tbRDA<-rda(as.formula(paste0("com.t","~",i)), data = env.t)
    
    pval<-permutest(tbRDA,10000)
    p.temp<-as.numeric(get.pval(pval))
    
    variance<-unname(tbRDA$CCA$eig/tbRDA$tot.chi*100)
    print(paste0(i," p:",round(p.temp,4)," var:",round(variance,1)))
  }
}

# import data
dat<-read.xlsx('Genus_Family_Order.xlsx', "Genus")
rownames(dat)<-dat[,1]
dat[,1]<-NULL
names(dat)<-gsub("X","",names(dat))

# import metadata
meta<-read.table(file = 'Metadata.csv', sep = ',', fileEncoding = 'UTF-8-BOM', stringsAsFactors = F, header = T)
rownames(meta)<-meta[,1]
meta[,1]<-NULL
meta[2:ncol(meta)] <- sapply(meta[2:ncol(meta)],as.numeric)

# exclude empty rows
dat<-dat[rowSums(dat)>0,]

# crosscheck with data
dat<-dat[,colnames(dat) %in% rownames(meta)]

# normalize counts 
dat.log <- dat + 1
dat.log <- log(dat.log)
dat.sc <- scale(dat.log, center = T, scale = F)
dat.sc <- as.data.frame(dat.sc)

# extra essens col
meta$veggy_true<-ifelse(meta$Meat_eater == 1 | meta$Fish_eater == 1, "0","1")
meta$veggy_true<-as.factor(meta$veggy_true)
meta$Meat_eater<-NULL
meta$Fish_eater<-NULL
meta$Vegan<-NULL
meta$Vegetarian<-NULL

# prepare for rda
env = meta
colnames(env)<-sub("-","",colnames(env))
colnames(env)<-sub("Tolerance","Tol",colnames(env))

# log transform meta
for(m in colnames(env[ ,!(colnames(env) %in% c("Sex", "Age","BMI","veggy_true","veggy_fish"))])){
  env[,m]<-env[,m]+(abs(min(env[,m],na.rm = T))+1)
  env[,m]<-log(env[,m],)
} 

# assign scaled data
com = as.data.frame(t(dat.sc))

# chose covariates you'd like to investigate
covars <- c("GCSF","MIP1a","MCP1","IL1RA","IL10","IL8","IL6","TNF","IP10")

# select them from data
env <- env[,covars]

# Create plot, pvals and explained variance
overview.plot2()
