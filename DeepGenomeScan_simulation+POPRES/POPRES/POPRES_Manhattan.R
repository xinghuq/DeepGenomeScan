library(vegan)
library(robust)
library(qvalue)
library(ade4)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(adegenet)
library(LEA)
library(pcadapt)
library(raster)
library(RColorBrewer)
library(splitstackshape)
library(sp)
library(gstat)
library(FactoMineR)
library(factoextra)
### read IMP for RD1 and RD2
read.excel <- function(header=TRUE,...) {
  read.table("clipboard",sep="\t",header=header,...)
}


IMP_KLFDAPC_RD1=read.csv("IMP_POPRES_mlpneuralnet_KLFDAPC5_RD1.csv")
IMP_KLFDAPC_RD2=read.csv("IMP_POPRES_mlpneuralnet_KLFDAPC5_RD2.csv")
### checking if Loci names are identical
identical(IMP_KLFDAPC_RD1$X,IMP_KLFDAPC_RD2$X)
POPRES_IMP_KLFDAPC_RD1_RD2=cbind(IMP_KLFDAPC_RD1,IMP_KLFDAPC_RD2$Overall)
save(POPRES_IMP_KLFDAPC_RD1_RD2,file = "POPRES_IMP_KLFDAPC_RD1_RD2.RData")
##For larger problems, the Orthogonalized Quadrant CorreRD2ion estimator is used.
DLqvaluesarsine<-function(DL_data,K)
{
  loadings<-DL_data# [,1:as.numeric(K)]
  normdat <- apply(loadings, 2, normalize)
  asindat=apply(normdat,2, function(x) {asin(sqrt(x))})
  resmaha <- covRob(asindat, distance = TRUE, na.action= na.omit, estim="donostah")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_DL<-qval$qvalues
  padj <- p.adjust(reschi2test,method="bonferroni")
  return(data.frame(p.values=reschi2test, q.values=q.values_DL,padj=padj,mahaD=resmaha))
}

###The covRob function selects a robust covariance estimator that is likely to provide a good estimate in a reasonable amount of time. Presently this selection is based on the problem size. The Donoho-Stahel estimator is used if there are less than 1000 observations and less than 10 variables or less than 5000 observations and less than 5 variables. If there are less than 50000 observations and less than 20 variables then the MCD is used. For larger problems, the Orthogonalized Quadrant Correlation estimator is used.
POPRES_KLFDAPC_pvalue=DLqvaluesarsine(POPRES_IMP_KLFDAPC_RD1_RD2[,-1],1)
POPRES_KLFDAPC_pvalue1=cbind(POPRES_IMP_KLFDAPC_RD1_RD2$X,POPRES_KLFDAPC_pvalue)
save(POPRES_KLFDAPC_pvalue1,file = "POPRES_KLFDAPC_pvalue.RData")

### read the loci position

POPTRS_loci_pos=read.table("GSK_08_24_01.EuroThinFinal.hg19.bim",header = FALSE)
POPTRS_loci_pos_anot=read.table("EuroThinFinal_dbSNP_anot_37.bim",header = FALSE)

##### align and identify the annotated 

match_pos1=POPTRS_loci_pos[match(POPRES_KLFDAPC_pvalue1$`POPRES_IMP_KLFDAPC_RD1_RD2$X`,POPTRS_loci_pos$V2),]

loci_pos_pvalue1=cbind(match_pos1,POPRES_KLFDAPC_pvalue1)
match_pos2=POPTRS_loci_pos_anot[match(loci_pos_pvalue1$V4,POPTRS_loci_pos_anot$V4),]

loci_pos_pvalue2=cbind(match_pos2,loci_pos_pvalue1)

save.image(file = "POPRES_DL_MLP_scan_KLFDAPC_pvalue_pos.RData")
## creat the ready to plot dataframe, BP, chr, p, snp names

POPRES_ManHattan_data=loci_pos_pvalue2[,c(1:2,4,14:16)]

colnames(POPRES_ManHattan_data)=c("CHR","SNP","BP","p.values","q.values","P")
##### making Manhattan plot 

POPRES_ManHattan_data_don=POPRES_ManHattan_data
##### update the latest annotated rs informatoion and genes



############################## annotating genes and SNP ids using chromsome position
## chr:position
posit1=read.excel()

##API key  2f0695fc0b019620fcc3217441765ecb7908


library(rentrez)
set_entrez_key("2f0695fc0b019620fcc3217441765ecb7908") ##please register your own account
Sys.getenv("ENTREZ_KEY")

entrez_db_searchable(db = "snp")

SNP_posit_gene=data.frame(matrix(ncol=6,nrow=220))
ano=list()
anodat=list()

for (i in seq_along(1: nrow(posit1))){
  res <- entrez_search(db = "snp", term = as.vector(posit1$pos)[i],version = c("2.0"),retmax = 3, sort="relevance",use_history = TRUE) 
  esums <- entrez_summary(db = "snp", id = min(as.integer(res$ids)))
  ano[[i]]=extract_from_esummary(esums, c("snp_id","chr","chrpos","genes","clinical_sort"))
  anodat[[i]]=do.call(cbind, ano[[i]])
  SNP_posit_gene[i,]=anodat[[i]]
}
library(data.table)
SNP_posit_gene1=do.call(rbind,ano)
SNP_posit_gene2=do.call(rbind,anodat)
SNP_posit_gene3=data.table::rbindlist(anodat,fill = TRUE)

write.csv(SNP_posit_gene1,file = "SNP_selected_annotated_rs_KLFDAPC1.csv")
write.csv(SNP_posit_gene3,file = "SNP_selected_annotated_rs_genes_KLFDAPC1.csv")


library(qqman)
Loci_under_selection=POPRES_ManHattan_data[(which(POPRES_ManHattan_data$P<1e-10)),]


write.csv(Loci_under_selection,file = "Loci_under_selection_KLFDAPC.csv")

save(Loci_under_selection,file = "Loci_under_selection_KLFDAPC.RData")



POPRES_ManHattan_data_don=POPRES_ManHattan_data
##### update the latest annotated rs informatoion and genes
Annotated_rs_genes=read.excel()   
POPRES_ManHattan_data_don$SNP=as.vector(POPRES_ManHattan_data_don$SNP)
Annotated_rs_genes$rsID=as.vector(Annotated_rs_genes$rsID)
POPRES_ManHattan_data_don[match(Annotated_rs_genes$BP,POPRES_ManHattan_data$BP),]$SNP=Annotated_rs_genes$rsID


save(POPRES_ManHattan_data_don,file = "POPRES_ManHattan_data_KLFDAPC.RData")

tiff(filename = "POPRES_manhattan_plot_KLFDAPC_1e_100_scale400_hightlight.tif",width = 10,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
manhattan(POPRES_ManHattan_data_don,ylim=c(0, 200),ylab=expression(-log[10](italic("q")-value)),annotatePval = 1e-10,highlight=Annotated_rs_genes$rsID,annotateTop = FALSE,suggestiveline = FALSE,genomewideline = -log10(1e-100))
dev.off()
pdf("POPRES_manhattan_plot_KLFDAPC_1e_100_scale400_hightlight.pdf",width = 10,height = 7,pointsize = 12)
manhattan(POPRES_ManHattan_data_don,ylim=c(0, 200),ylab=expression(-log[10](italic("q")-value)),annotatePval = 1e-10,highlight=Annotated_rs_genes$rsID,annotateTop = FALSE,suggestiveline = FALSE,genomewideline = -log10(1e-100))
dev.off()

tiff(filename = "POPRES_manhattan_KLFDAPC_plot_1e_100_scale300_ylab_hightlight.tif",width = 10,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
manhattan(POPRES_ManHattan_data_don,ylim=c(0, 200),ylab=expression(-log[10](italic("q")-value)),annotatePval = 1e-10,highlight=Annotated_rs_genes$rsID,annotateTop = FALSE,suggestiveline = FALSE,genomewideline = -log10(1e-100))
dev.off()
pdf("POPRES_manhattan_plot_KLFDAPC_1e_100_scale300_ylab_hightlight.pdf",width = 10,height = 7,pointsize = 12)
manhattan(POPRES_ManHattan_data_don,ylim=c(0, 200),ylab=expression(-log[10](italic("q")-value)),annotatePval = 1e-10,highlight=Annotated_rs_genes$rsID,annotateTop = FALSE,suggestiveline = FALSE,genomewideline = -log10(1e-100))
dev.off()

save.image(file = "POPRES_DL_MLP_scan_KLFDAPC_pvalue_Manhanton_plot.RData")

save(POPRES_ManHattan_data_don,file = "POPRES_ManHattan_data_don_KLFDAPC_R1_R2.RData")

