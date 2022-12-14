# This script is to generate table of correlation analysis between expression level of 410 RBP and 377 miRNA by each site primary
# 2022/12/06 made

# make new directory
setwd("C:/Rdata")
dir.create("20221206_correlaiton_analysis_between_RBP_and_miRNA_in_each_site")

# import result of correlation analysis
# this result is located at "https://github.com/Ryosuke-Hirota/20221122_ROC_curve_of_pvalue_after_cutoff"
setwd("C:/Rdata/20221108_correlaiton_between_RBP_level_and_miRNA_level_with_excluding_outlier")
cor.result <-read.table("summary_of_correlaiton_between_RBP_level_and_miRNA_level_with_excluding_outlier_and_considering_no_expression.txt",sep="\t",header = T,stringsAsFactors = F)
outliers <-read.table("table_about_outliers_in_correlaiton_analysis_between_RBP_level_and_miRNA_level.txt",sep="\t",header = T,stringsAsFactors = F)
colnames(outliers)[6] <-"miRNA_outlier_large"
cor.out.df <-merge(cor.result,outliers,by=c("RBP","miRNA"))

# import table about expression level of RBP
# this table is located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\CCLE_Data"
setwd("C:/Rdata/CCLE_data")
RBP.df <-read.table("RNAseq_of_RBP.txt",sep="\t",header = T,stringsAsFactors = F)
RBP.df <-RBP.df[,c(2:411,1,412)]

# get name of site primary
site <-sort(unique(RBP.df[,412]))

# import table about expression level of miRNAs
# this table is located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\CCLE_Data"
miRNA.df <-read.table("CCLE_miRNAseq.txt",sep="\t",header = T,stringsAsFactors = F,check.names = F)
miRNA.df <-miRNA.df[,c(2:378,1,379)]

# set function for calculating minimum value of expression level
find.min <-function(x,y){
  for (i in 1:y) {
    m <-x[x[,i]!=0,i]
    mv <-min(m)
    if(i==1){
      mvs <-mv
    }else{
      mvs <-append(mvs,mv)
    }}
  return(min(mvs))
}

# investigate minimum value of each expression level
r.min <-find.min(RBP.df,410)
m.min <-find.min(miRNA.df,377)

# make empty lists to put in each expression level in each site
site.RBP.list <-as.list(NULL)
site.miRNA.list <-as.list(NULL)

# list expression level of RBP or miRNA by each site
for (s in 1:length(site)) {
  site.RBP.df <-RBP.df[RBP.df[,412]==site[s],]
  site.miRNA.df <-miRNA.df[miRNA.df[,379]==site[s],]
  site.RBP.list[[s]] <-site.RBP.df
  site.miRNA.list[[s]] <-site.miRNA.df
  }

sm <-as.data.frame(matrix(nrow = 410*377,ncol = 7))
colnames(sm) <-c("RBP","miRNA","r","p.value","number_of_cell_line","number_of_cell_line_before_excluding","number_of_all_cell_line")

setwd("C:/Rdata/20221206_correlaiton_analysis_between_RBP_and_miRNA_in_each_site")

# generate summary
for (s in 1:length(site)) {
  for (i in 1:410) {
    for (k in 1:377) {
    
    # extract expression levels of RBP and miRNA
    ex.site.RBP <-as.data.frame(site.RBP.list[[s]])
    RBP <-ex.site.RBP[,c(i,411,412)]
    ex.site.miRNA <-as.data.frame(site.miRNA.list[[s]])
    miRNA <-ex.site.miRNA[,c(k,378,379)]
    RBP.miRNA <-merge(RBP,miRNA,by=c("CELL","Site_Primary"))
    RBP.miRNA <-RBP.miRNA[,c(3,4,1,2)]
    
    # count number of row
    site.all.cell <-nrow(RBP.miRNA)
    
    #
    combination <-cor.out.df[cor.out.df[,1]==colnames(RBP.miRNA)[1]&cor.out.df[,2]==colnames(RBP.miRNA)[2],]
    zero.cell <-combination[,6]
    
    # if number of cell line without RBP expression is higher than 100 in result of correlation analysis, add minimum RBP expression level
    # if number of cell line without RBP expression is lower than 100 in result of correlation analysis, exclude cell line without RBP expression
    if(zero.cell>100){
      RBP.miRNA[,1] <-ifelse(RBP.miRNA[,1]==0,r.min,RBP.miRNA[,1])
    }else{
      RBP.miRNA <-RBP.miRNA[RBP.miRNA[,1]!=0,]
    }
    
    # log2 each expression level 
    RBP.miRNA[,1] <-log2(RBP.miRNA[,1])
    RBP.miRNA[,2] <-log2(RBP.miRNA[,2])

    
    # exclude outlier of RBP expression level (outlier was determined in correlation analysis between 410RBPs and 377 miRNAs)
    RBP.miRNA <-RBP.miRNA[RBP.miRNA[,1]>=combination[,7],]
    RBP.miRNA <-RBP.miRNA[RBP.miRNA[,1]<=combination[,8],]
    
    # exclude outlier of miRNA expression level (outlier was determined in correlation analysis between 410RBPs and 377 miRNAs)
    RBP.miRNA <-RBP.miRNA[RBP.miRNA[,2]>=combination[,9],]
    RBP.miRNA <-RBP.miRNA[RBP.miRNA[,2]<=combination[,10],]
    
    # calculate correlation coefficient
    r <-try(cor.test(RBP.miRNA[,1],RBP.miRNA[,2],method = "pearson"),silent = T)
    
    if(class(r)!="try-error"){
    # write summary
    sm[(i-1)*377+k,1] <-colnames(RBP.miRNA)[1]
    sm[(i-1)*377+k,2] <-colnames(RBP.miRNA)[2]
    sm[(i-1)*377+k,3] <-signif(r$estimate,3)
    sm[(i-1)*377+k,4] <-signif(r$p.value,3)
    sm[(i-1)*377+k,5] <-nrow(RBP.miRNA)
    sm[(i-1)*377+k,6] <-site.all.cell
    sm[(i-1)*377+k,7] <-combination[,5]
    }else{
      sm[(i-1)*377+k,1] <-colnames(RBP.miRNA)[1]
      sm[(i-1)*377+k,2] <-colnames(RBP.miRNA)[2]
      sm[(i-1)*377+k,3] <-NA
      sm[(i-1)*377+k,4] <-NA
      sm[(i-1)*377+k,5] <-nrow(RBP.miRNA)
      sm[(i-1)*377+k,6] <-site.all.cell
      sm[(i-1)*377+k,7] <-combination[,5]
      }
    
    if(i==410&k==377){
      # output summary
      write.table(sm,paste0(site[s],"_summary_of_correlaiton_between_RBP_level_and_miRNA_level_with_excluding_outlier_and_considering_no_expression.txt"),
                            sep="\t",row.names = F,quote = F)
    }else{}
    }}}

# list tables of each site
sum.list <-list.files(path="C:/Rdata/20221206_correlaiton_analysis_between_RBP_and_miRNA_in_each_site",pattern = ".txt")

# make table summarized tables of each site
for (l in 1:length(sum.list)) {
  df <-read.table(sum.list[l],sep="\t",header = T,stringsAsFactors = F)
  df[,1] <-paste0(df[,1],"_vs_",df[,2])
  df <-df[,c(1,3)]
  colnames(df) <-c("combination",site[l])
  if(l==1){
    cor.table <-df
  }else{
    cor.table <-merge(cor.table,df,by="combination")
  }
}

# output summary
write.table(cor.table,"table_of_correlation_coefficient_by_each_site.txt",sep="\t",row.names = F,quote = F)
