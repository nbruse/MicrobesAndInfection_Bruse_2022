# Analysis Microbiome Bruse et al.
# Author: NB 

library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(xlsx)
library(tibble)
library(gtools)

# corr alpha cytokines
setwd("")

clean<-function(x){
  x<-gsub("_D0_T0", "", x)
  x<-gsub("_D0", "", x)
  x<-gsub("Max.delta", "", x)
  x<-gsub("PEAK_", "", x)
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
clean.m<-function(x){
  x<-gsub("\\[", "", x)
  x<-gsub("]", "", x)
  x<-gsub("_", " ", x)
  x<-gsub("R-7", "R7", x)
  x<-gsub('X.',"", x)
  
  return(x)
}

# read metadata
meta<-read.csv("Metadata.csv", stringsAsFactors = F)
rownames(meta)<-meta[,1]
meta[,1]<-NULL
meta$Vegetarian<-NULL
meta$Vegan<-NULL
meta$Fish_eater<-NULL
meta$Meat_eater<-NULL
meta[2:ncol(meta)] <- sapply(meta[2:ncol(meta)],as.numeric)

# import data
dat<-read.xlsx("Genus_Family_Order.xlsx", "Genus")
rownames(dat)<-dat[,1]
dat[,1]<-NULL
names(dat)<-gsub("X","", names((dat)))
dat<-dat[,colnames(dat) %in% rownames(meta)]
dat$median<-apply(dat, 1, median)
dat2 <-dat %>% 
  rownames_to_column('biota') %>% 
  arrange(desc(median)) %>% 
  column_to_rownames('biota') %>%
  filter(rownames(.) != "UCG-002") %>%
  top_n(25) %>% 
  select(-median)

#### Pairwise correlation ####
# select ev
ev <- meta[,22:61]
dat2<-data.frame(t(dat2))

# make rhos en pvals
l_p<-list()
l_r<-list()
for(i in colnames(ev)){
  pvals<-c()
  rhos<-c()
  
  for(h in colnames(dat2)){
    
    test<-cor.test(dat2[,h], ev[,i], method = "spearman", exact = F, correct = F)
    
    p<-test$p.value
    rho<-unname(test$estimate)
    
    pvals<-c(pvals,p)
    rhos<-c(rhos,rho)
  }
  l_p[[i]]<-pvals
  l_r[[i]]<-rhos
}

out_p<-do.call(cbind.data.frame,l_p)
out_p<-as.data.frame(t(apply(out_p, 1, p.adjust, method = "fdr")))
rownames(out_p)<-colnames(dat2)

# make rho df
out_r<-do.call(cbind.data.frame,l_r)
rownames(out_r)<-colnames(dat2)

# clean up rho df
out_r[is.na(out_p)] <- NA

# prep pvals to stars
out_p<-apply(out_p, 1, stars.pval) %>% as.data.frame
out_p<-as.data.frame(lapply(out_p, function(y) gsub("\\.", "", y)))
rownames(out_p)<-names(l_p)

# clean colnames 
colnames(out_r)<-clean(colnames(out_r))
rownames(out_r)<-clean.m(rownames(out_r))
colnames(out_p)<-clean.m(colnames(out_p))

# make heatmap 
out_r %>% 
  t() %>% 
  as.data.frame() %>%
  pheatmap(cluster_rows = F, cluster_cols = F, na_col = "white", color = rev(brewer.pal(n = 11, "RdBu")),
           breaks = seq(-.5,.5, length = 11), border_color = 'black',
           cellwidth = 20, cellheight = 13,
           display_numbers = out_p,fontsize_number = 17, number_color = 'white', fontsize = 15,angle_col = 315)

#### Correlation with alpha ####
# import alpha index file
alpha<-read.table(file = "alpha-diversity_SHANNON.tsv", sep ="\t", header = T,row.names = 1)
alpha<-alpha[rownames(alpha) %in% rownames(meta),,drop=F]

# make Rhos en pvals
pvals<-c()
rhos<-c()
names_rho<-c()
for(i in colnames(meta)){
  if(is.numeric(meta[,i])){
    test<-cor.test(alpha$shannon, meta[rownames(meta) %in% rownames(alpha),i], method = "spearman", exact = F, correct = F)
    p<-test$p.value
    rho<-unname(test$estimate)
    names_rho<-c(names_rho,i)
    pvals<-c(pvals,p)
    rhos<-c(rhos,rho)
  }
}

# make rho df
rhos[is.na(pvals)]<-NA
names(rhos)<-names_rho
names(pvals)<-names_rho

# prepare pval df for heatmap
pvals2<- pvals %>% t() %>% data.frame() %>%
  t() %>% data.frame() %>%
  tibble::rownames_to_column("names") %>%
  dplyr::rename('value'='.') %>%
  mutate(value=stars.pval(value)) %>%
  mutate(value=gsub("\\.", "", value))

# prepare rho df for heatmap
heat<-rhos %>% t() %>% data.frame() %>%
  t() %>% data.frame() %>%
  tibble::rownames_to_column("names") %>%
  left_join(pvals2, by = "names") %>%
  tibble::column_to_rownames('names') %>%
  dplyr::rename('Shannon index'='.')

values<-heat[,1,drop=F]
ps<-heat[,2,drop=F]
ps$value[is.na(ps$value)]<-''

# clean names
rownames(values)<-clean(rownames(values))

# heatmap
pheatmap(values,cluster_rows = F,cluster_cols = F, na_col = "white", fontsize = 20, color = rev(brewer.pal(n = 11, "RdBu")),
         display_numbers = ps, fontsize_number = 25, cellheight = 22,number_color = 'white',angle_col = 0)
