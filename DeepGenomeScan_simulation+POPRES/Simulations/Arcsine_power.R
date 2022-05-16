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


#### These are scripts for analysis of simulations using three approaches, PCAdapt, RDA, and DeepGenomeScan
### read the importance files produced from deepgenomescan
impdir=list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
impfile <- list.files(impdir, pattern="*_mlp_olden2_importance_env", full.names = TRUE)
## 100 simulations and 10 envir variables, 

impSim_all=lapply(impfile, read.csv) ### including all 100 simulations and 10 variables each simuilations, 
impfile[91:100]## check, each 10 file is a simulation, 18
### x is the loci ID and Overall is the overall importance, 1000 list: including 100 simulations and 10 envir
impfile[71:80]#
### combine every factor (envir) into one file, meaning combining each single simulation into one file

## check impSim_all[[11:20]]### This is sim1, check the list file
imp_sim_list=as.data.frame(do.call(cbind, impSim_all))
### remove redundant loci ID
Xrep<- !names(imp_sim_list) %in% c("X")
DLdata <- imp_sim_list[,Xrep, drop = FALSE]
colnames(DLdata) ### now 1:10 is the first simulation
### spilt into 100 simulations each by 10
save.image(file = "DL_Scan2.RData")

### separate each simulation

DLsim1=DLdata[71:80] ## This is the simulation 1, check
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
chunk(1:1000,100)
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
chunk2(1:1000,100)
samplesim=chunk2(1:1000,100)

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

Loci<-rep("Neutral", 1000)
Loci[c(1,11,21,31,41,51,61,71,81,91)]<-"QTL1"
Loci[c(101,111,121,131,141,151,161,171,181,191)]<-"QTL2"
Loci[c(201,211,221,231,241,251,261,271,281,291)]<-"QTL3"
Selected_Loci<-Loci[-which(Loci=="Neutral")]


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
DLsim1=apply(DLdata[91:100],2,normalize) #18
# DLsim1=apply(DLdata[71:80],2,normalize)# 16
Simqvaluear=DLqvaluesarsine(DLsim1,10)



## Manhattan plot

ggplot() +
  geom_point(aes(x=which(Loci=="Neutral"), y=-log10(Simqvaluear[-which(Loci!="Neutral"),1])), col = "gray83") +
  geom_point(aes(x=which(Loci!="Neutral"), y=-log10(Simqvaluear[-which(Loci=="Neutral"),1]), colour = Selected_Loci)) +
  xlab("SNPs") + ylab("DNN -log10(p-value)") +ylim(c(0,100))+theme_bw()



#### This part is a Genome Scan using RDA


#### Function of RDA-based genome scan ####

rdadapt<-function(rda,K)
{
  loadings<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(loadings, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

##############
#### Data ####

# Genet
geno<-t(read.table("geno.pcadapt"))

geno<-read.csv("sim1.csv")[,-c(1:14)]
sp <- apply(geno, 2, mean) ### column mean--- loci mean number of ancesteral alleles
geno <-  geno[ , sp > 0.05 & sp < 1.95]

# Env
env<-read.csv("sim1.csv")
env<-as.data.frame(apply(env, 2, scale))


#############################
#### RDA on simulation 1 ####

RDA_tot <- rda(geno ~ env$envir1 + env$envir2 + env$envir3 + env$envir4 + env$envir5 + env$envir6 + env$envir7 + env$envir8 + env$envir9 + env$envir10,  env)
res_rdadapt<-rdadapt(RDA_tot, 10)

## Plot of SNPs and variables in the RDA space
coord_snps <- RDA_tot$CCA$v[, 1:5]
coord_env <- RDA_tot$CCA$biplot[, 1:5]
row.names(coord_env) <- colnames(env[,2:11])

Loci<-rep("Neutral", 1000)
Loci[c(1,11,21,31,41,51,61,71,81,91)]<-"QTL1"
Loci[c(101,111,121,131,141,151,161,171,181,191)]<-"QTL2"
Loci[c(201,211,221,231,241,251,261,271,281,291)]<-"QTL3"
Selected_Loci<-Loci[-which(Loci=="Neutral")]

p1 <- ggplot() +
  geom_point(aes(x=coord_snps[-which(Loci!="Neutral"),1], y=coord_snps[-which(Loci!="Neutral"),2]), col = "gray86") +
  geom_point(aes(x=coord_snps[which(Loci!="Neutral"),1], y=coord_snps[which(Loci!="Neutral"),2], colour = Selected_Loci), size = 2.5) +
  geom_segment(aes(xend=coord_env[,1]/10, yend=coord_env[,2]/10, x=0, y=0), colour="black", size=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(aes(x=1.2*coord_env[,1]/10, y=1.2*coord_env[,2]/10, label = row.names(coord_env)), size=3) +
  xlab("RDA 1") + ylab("RDA 2") +
  theme_bw() +
  theme(legend.position="none")

p2 <- ggplot() +
  geom_point(aes(x=coord_snps[-which(Loci!="Neutral"),1], y=coord_snps[-which(Loci!="Neutral"),3]), col = "gray86") +
  geom_point(aes(x=coord_snps[which(Loci!="Neutral"),1], y=coord_snps[which(Loci!="Neutral"),3], colour = Selected_Loci), size = 2.5) +
  geom_segment(aes(xend=coord_env[,1]/10, yend=coord_env[,3]/10, x=0, y=0), colour="black", size=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(aes(x=1.2*coord_env[,1]/10, y=1.2*coord_env[,3]/10, label = row.names(coord_env)), size=3) +
  xlab("RDA 1") + ylab("RDA 3") +
  theme_bw() +
  theme(legend.position="none")

p3 <- ggplot() +
  geom_point(aes(x=coord_snps[-which(Loci!="Neutral"),1], y=coord_snps[-which(Loci!="Neutral"),4]), col = "gray86") +
  geom_point(aes(x=coord_snps[which(Loci!="Neutral"),1], y=coord_snps[which(Loci!="Neutral"),4], colour = Selected_Loci), size = 2.5) +
  geom_segment(aes(xend=coord_env[,1]/10, yend=coord_env[,4]/10, x=0, y=0), colour="black", size=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(aes(x=1.2*coord_env[,1]/10, y=1.2*coord_env[,4]/10, label = row.names(coord_env)), size=3) +
  xlab("RDA 1") + ylab("RDA 4") +
  theme_bw()

plot_grid(p1, p2, p3, labels = "AUTO", align = "h", rel_widths = c(.6, .6, 1), ncol=3)


## Manhattan plot
ggplot() +
  geom_point(aes(x=which(Loci=="Neutral"), y=-log10(res_rdadapt[-which(Loci!="Neutral"),1])), col = "gray83") +
  geom_point(aes(x=which(Loci!="Neutral"), y=-log10(res_rdadapt[-which(Loci=="Neutral"),1]), colour = Selected_Loci)) +
  xlab("SNPs") + ylab("RDADAPT with condition)") +
  theme_bw()


#################################
#### PCADAPT on simulation 1 ####

setwd("./sim1/")
data<-read.pcadapt("geno.pcadapt", type = "pcadapt")

pcadapt <- pcadapt(data, K=10, method="mahalanobis",  min.maf = 0.05, ploidy = 2)


##########################################
#### Comparison PCADAPT - RDA on sim1 ####

## Screeplots
p1 <- ggplot() +
  geom_line(aes(x=c(1:10), y=as.vector(pcadapt$singular.values)), linetype="dotted", size = 1.5, colour="black") +
  geom_point(aes(x=c(1:10), y=as.vector(pcadapt$singular.values)), size = 4, colour="black", shape=17) +
  scale_x_discrete(name = "Ordination axes", limits=c(1:10)) +
  ylab("Inertia") +
  ggtitle("PCA") +
  theme_bw()
p2 <- ggplot() +
  geom_line(aes(x=c(1:10), y=as.vector(RDA_tot$CCA$eig)), linetype="dotted", size = 1.5, colour="grey") +
  geom_point(aes(x=c(1:10), y=as.vector(RDA_tot$CCA$eig)), size = 4, colour="grey", shape=19) +
  scale_x_discrete(name = "Ordination axes", limits=c(1:10)) +
  ylab("Inertia") +
  ggtitle(label = "RDA") +
  theme_bw()
ggarrange(p1, p2, ncol=2, nrow = 1)

## Manhattan plot PCA & RDA with 10 axes
pcadapt <- pcadapt(data, K=10, method="mahalanobis",  min.maf = 0.05, ploidy = 2)

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}


#donostah is good if variables is less than 1000

Loci<-rep("Neutral", 1000)
Loci[c(1,11,21,31,41,51,61,71,81,91)]<-"QT1"
Loci[c(101,111,121,131,141,151,161,171,181,191)]<-"QT2"
Loci[c(201,211,221,231,241,251,261,271,281,291)]<-"QT3"
### p values

p.values<-data.frame(index = c(1:1000), pcadapt = -log10(pcadapt$pvalues), rdadapt = -log10(res_rdadapt[,1]), DL_mlp=-log10(Simqvaluear$p.values),Loci=Loci)

Selected_Loci<-p.values$Loci[-which(p.values$Loci=="Neutral")]

## QQ plot DL_SCAN
tiff(filename = "QQ_plot_DL_SCAN_Arc-sin_small_letters.tif",width = 7.5,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
qqplot(rexp(length(Simqvaluear[,1]), rate = log(10)),
       -log10(Simqvaluear[,1]), xlab = "Expected quantile",ylab = "-log10(p-value)",
       pch = 19, cex = .4)
abline(0,1)
dev.off()
pdf("QQ_plot_DL_SCAN_Arc-sin_small_letters.pdf",width = 7.5,height = 7,pointsize = 12)
qqplot(rexp(length(Simqvaluear[,1]), rate = log(10)),
       -log10(Simqvaluear[,1]), xlab = "Expected quantile",ylab = "-log10(p-value)",
       pch = 19, cex = .4)
abline(0,1)
dev.off()

### QQ plot pcadapt
tiff(filename = "QQ_plot_PCAdapt_small_letters.tif",width = 7.5,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
qqplot(rexp(length(pcadapt$pvalues), rate = log(10)),
       -log10(pcadapt$pvalues), xlab = "Expected quantile",ylab = "-log10(p-value)",
       pch = 19, cex = .4)
abline(0,1)
dev.off()
pdf("QQ_plot_PCAdapt_small_letters.pdf",width = 7.5,height = 7,pointsize = 12)
qqplot(rexp(length(pcadapt$pvalues), rate = log(10)),
       -log10(pcadapt$pvalues), xlab = "Expected quantile",ylab = "-log10(p-value)",
       pch = 19, cex = .4)
abline(0,1)
dev.off()

### QQ plot RDA
tiff(filename = "QQ_plot_RDA_small_letters.tif",width = 7.5,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
qqplot(rexp(length(res_rdadapt$p.values), rate = log(10)),
       -log10(res_rdadapt$p.values), xlab = "Expected quantile",ylab = "-log10(p-value)",
       pch = 19, cex = .4)
abline(0,1)
dev.off()
pdf("QQ_plot_RDA_small_letters.pdf",width = 7.5,height = 7,pointsize = 12)
qqplot(rexp(length(res_rdadapt$p.values), rate = log(10)),
       -log10(res_rdadapt$p.values), xlab = "Expected quantile",ylab = "-log10(p-value)",
       pch = 19, cex = .4)
abline(0,1)
dev.off()
#### Fig. S3


S3p1 <- ggplot() +
  geom_point(aes(x=p.values$index[-which(p.values$Loci!="Neutral")], y=p.values$pcadapt[-which(p.values$Loci!="Neutral")]), col = "gray83") +
  geom_point(aes(x=as.vector(p.values$index[-which(p.values$Loci=="Neutral")]), y=as.vector(p.values$pcadapt[-which(p.values$Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNP (with maf > 0.05)") + ylab("PCADAPT -log10(p-values)") +
  ylim(0,10) +
  theme_bw()
S3p2 <- ggplot() +
  geom_point(aes(x=p.values$index[-which(p.values$Loci!="Neutral")], y=p.values$rdadapt[-which(p.values$Loci!="Neutral")]), col = "gray83") +
  geom_point(aes(x=as.vector(p.values$index[-which(p.values$Loci=="Neutral")]), y=as.vector(p.values$rdadapt[-which(p.values$Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNP (with maf > 0.05)") + ylab("RDA -log10(p-values)") +
  ylim(0,10) +
  theme_bw()
S3p3 <- ggplot() +
  geom_point(aes(x=p.values$index[-which(p.values$Loci!="Neutral")], y=p.values$DL_mlp[-which(p.values$Loci!="Neutral")]), col = "gray83") +
  geom_point(aes(x=as.vector(p.values$index[-which(p.values$Loci=="Neutral")]), y=as.vector(p.values$DL_mlp[-which(p.values$Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNP (with maf > 0.05)") + ylab("DNN -log10(p-values)") +
  ylim(0,100) +
  theme_bw()

ggarrange(S3p1, S3p2, S3p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE, legend="right")

tiff(filename = "PCA_RDA_DL_MLP_Scan_Sim1_qvalue_all1_small_letters20220419.tif",width = 21,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
ggarrange(S3p1, S3p2, S3p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE, font.label = list(size=18),  legend="right")
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_Sim1_qvalue_all1_small_letters20220419.pdf",width = 21,height = 7,pointsize = 12)
ggarrange(S3p1, S3p2, S3p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE, font.label = list(size=18),  legend="right")
dev.off()

tiff(filename = "PCA_RDA_DL_MLP_Scan_Sim1_14_7_qvalue_all1_small_letters20220419.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
ggarrange(S3p1, S3p2, S3p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE,font.label = list(size=18),   legend="right")
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_Sim1_14_7_qvalue_all1_small_letters20220419.pdf",width = 14,height = 7,pointsize = 12)
ggarrange(S3p1, S3p2, S3p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE,font.label = list(size=18),   legend="right")
dev.off()



p1 <- ggplot() +
  geom_point(aes(x=p.values$index[-which(p.values$Loci!="Neutral")], y=p.values$pcadapt[-which(p.values$Loci!="Neutral")]), col = "gray83") +
  geom_point(aes(x=as.vector(p.values$index[-which(p.values$Loci=="Neutral")]), y=as.vector(p.values$pcadapt[-which(p.values$Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNP (with maf > 0.05)") + ylab("PCADAPT -log10(q-values)") +
  ylim(0,100) +
  theme_bw()
p2 <- ggplot() +
  geom_point(aes(x=p.values$index[-which(p.values$Loci!="Neutral")], y=p.values$rdadapt[-which(p.values$Loci!="Neutral")]), col = "gray83") +
  geom_point(aes(x=as.vector(p.values$index[-which(p.values$Loci=="Neutral")]), y=as.vector(p.values$rdadapt[-which(p.values$Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNP (with maf > 0.05)") + ylab("RDA_Scan -log10(q-values)") +
  ylim(0,100) +
  theme_bw()
p3 <- ggplot() +
  geom_point(aes(x=p.values$index[-which(p.values$Loci!="Neutral")], y=p.values$DL_mlp[-which(p.values$Loci!="Neutral")]), col = "gray83") +
  geom_point(aes(x=as.vector(p.values$index[-which(p.values$Loci=="Neutral")]), y=as.vector(p.values$DL_mlp[-which(p.values$Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNP (with maf > 0.05)") + ylab("DNN -log10(q-values)") +
  ylim(0,100) +
  theme_bw()

#### 05-09-2020 update change front letters as small letters


tiff(filename = "PCA_RDA_DL_MLP_Scan_Sim1_log_qvalues1_small_letters20220419_.tif",width = 21,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
ggarrange(p1, p2, p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE,font.label = list(size=18), legend="right")
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_Sim1_log_qvalues1_small_letters20220409_.pdf",width = 21,height = 7,pointsize = 12)
ggarrange(p1, p2, p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE, font.label = list(size=18), legend="right")
dev.off()

tiff(filename = "PCA_RDA_DL_MLP_Scan_Sim1_14_7_log_qvalues1_small_letters20220409_.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
ggarrange(p1, p2, p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE,font.label = list(size=18),  legend="right")
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_Sim1_14_7_log_qvalues1_small_letters20220409_.pdf",width = 14,height = 7,pointsize = 12)
ggarrange(p1, p2, p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE, font.label = list(size=18), legend="right")
dev.off()
########


#######################################################################################################
#### Comparison DL_mlp PCADAPT & RDA over the 100 simulations ####
#### DL_mlp
### split the data into 100 simulations
DLdatasim=lapply(samplesim,function(x) {DLdata[,x]})
p.values_DLmlp_simusarcsin <- lapply(DLdatasim, function(x) {DLqvaluesarsine(x, 10)[,1]})


### adj.pvalue is another form of q.value

#RDA
fich_simu<-list.files(recursive = T, path="./", pattern="geno.pcadapt", full.names = TRUE)
genos<-lapply(fich_simu, function(x) t(read.table(x)))
env<-read.csv("./sim1/sim1.csv")#[,1:10]

RDA_res <- lapply(genos, function(x) {rda(x ~ env$envir1 + env$envir2 + env$envir3 + env$envir4 + env$envir5 + env$envir6 + env$envir7 + env$envir8 + env$envir9 + env$envir10,  env)})
### the original article k =5
p.values_rdadapt_simus <- lapply(RDA_res, function(x) {rdadapt(x, 10)[,1]})
# PCADAPT
genos_pcadapt<-lapply(fich_simu, function(x) read.pcadapt(x, type = "pcadapt"))
pcadaptsimall <- lapply(genos_pcadapt, function(x) pcadapt(x, K=10, method="mahalanobis",  min.maf = 0.05, ploidy = 2))
p.values_pcadapt_simus <- lapply(pcadaptsimall, function(x) {x$pvalues})

save.image(file = "DL_Scan_PCA_RDA_DL_mlp_100sim.RData")
## Power = number of true positive found by the analysis for each QTL
Loci<-rep("Neutral", 1000)
Loci[c(1,11,21,31,41,51,61,71,81,91)]<-"QT1"
Loci[c(101,111,121,131,141,151,161,171,181,191)]<-"QT2"
Loci[c(201,211,221,231,241,251,261,271,281,291)]<-"QT3"
### first set the cutoff as 0.01, number of loci detected by a method per sim

#### We set p value for RDA and PCAdapt as 1e-8 to contrl FDR
Tpos_pcadapt_QTL1_e8<-unlist(lapply(p.values_pcadapt_simus, function(x) length(which(x < 1e-8 & Loci=="QT1"))))
Tpos_pcadapt_QTL2_e8<-unlist(lapply(p.values_pcadapt_simus, function(x) length(which(x < 1e-8 & Loci=="QT2"))))
Tpos_pcadapt_QTL3_e8<-unlist(lapply(p.values_pcadapt_simus, function(x) length(which(x < 1e-8 & Loci=="QT3"))))

Tpos_rdadapt_QTL1_e8<-unlist(lapply(p.values_rdadapt_simus, function(x) length(which(x < 1e-8 & Loci=="QT1"))))
Tpos_rdadapt_QTL2_e8<-unlist(lapply(p.values_rdadapt_simus, function(x) length(which(x < 1e-8 & Loci=="QT2"))))
Tpos_rdadapt_QTL3_e8<-unlist(lapply(p.values_rdadapt_simus, function(x) length(which(x < 1e-8 & Loci=="QT3"))))

Tpos_DLmlp_QTL1_e8<-unlist(lapply(p.values_DLmlp_simusarcsin, function(x) length(which(x < 1e-8  & Loci=="QT1"))))
Tpos_DLmlp_QTL2_e8<-unlist(lapply(p.values_DLmlp_simusarcsin, function(x) length(which(x < 1e-8 & Loci=="QT2"))))
Tpos_DLmlp_QTL3_e8<-unlist(lapply(p.values_DLmlp_simusarcsin, function(x) length(which(x < 1e-8  & Loci=="QT3"))))

## adju pvalue, do these later to compare the p.values and then choose one form


## False postives 
Fpos_pcadapt_e8<-unlist(lapply(p.values_pcadapt_simus, function(x) length(which(x < 1e-8 & Loci=="Neutral"))))
Fpos_rdadapt_e8<-unlist(lapply(p.values_rdadapt_simus, function(x) length(which(x < 1e-8 & Loci=="Neutral"))))

Fpos_DLmlp_e8<-unlist(lapply(p.values_DLmlp_simusarcsin, function(x) length(which(x <1e-8 & Loci=="Neutral"))))


### t test to see the difference
t.test(Tpos_DLmlp_QTL1_e8/10,Tpos_rdadapt_QTL1_e8/10)
t.test(Tpos_DLmlp_QTL1_e8/10,Tpos_pcadapt_QTL1_e8/10)
t.test(Tpos_rdadapt_QTL1_e8/10,Tpos_pcadapt_QTL1_e8/10)


t.test(Tpos_DLmlp_QTL2_e8/10,Tpos_rdadapt_QTL2_e8/10)
t.test(Tpos_DLmlp_QTL2_e8/10,Tpos_pcadapt_QTL2_e8/10)
t.test(Tpos_rdadapt_QTL2_e8/10,Tpos_pcadapt_QTL2_e8/10)

t.test(Tpos_DLmlp_QTL3_e8/10,Tpos_rdadapt_QTL3_e8/10)
t.test(Tpos_DLmlp_QTL3_e8/10,Tpos_pcadapt_QTL3_e8/10)
t.test(Tpos_rdadapt_QTL3_e8/10,Tpos_pcadapt_QTL3_e8/10)


Allpositivespcadapt_e8=unlist(lapply(p.values_pcadapt_simus, function(x) length(which(x < 1e-8))))
Allpositivesrda_e8=unlist(lapply(p.values_rdadapt_simus, function(x) length(which(x < 1e-8))))
AllpositivesDLmlp_e8=unlist(lapply(p.values_DLmlp_simusarcsin, function(x) length(which(x < 1e-8))))


### True Discovery rate
Tpos_e8<-data.frame(Tpos=c(Tpos_pcadapt_QTL1_e8/10, Tpos_pcadapt_QTL2_e8/10, Tpos_pcadapt_QTL3_e8/10, Tpos_rdadapt_QTL1_e8/10, Tpos_rdadapt_QTL2_e8/10, Tpos_rdadapt_QTL3_e8/10,Tpos_DLmlp_QTL1_e8/10, Tpos_DLmlp_QTL2_e8/10, Tpos_DLmlp_QTL3_e8/10), Method=c(rep("PCA", 300), rep("RDA", 300),rep("DNN", 300)), Loci=c(rep("QT1", 100), rep("QT2", 100), rep("QT3", 100), rep("QT1", 100), rep("QT2", 100), rep("QT3", 100),rep("QT1", 100), rep("QT2", 100), rep("QT3", 100)))
Tpos_meansd_e8<-cbind(aggregate(Tpos_e8$Tpos, by=list(Tpos_e8$Method,Tpos_e8$Loci), FUN = "mean"), sd=aggregate(Tpos_e8$Tpos, by=list(Tpos_e8$Method,Tpos_e8$Loci), FUN = "sd")$x)
colnames(Tpos_meansd_e8)[1:3]<-c("Method", "Loci", "Rate")




## Visualisation
Tposplot <- ggplot(Tpos_meansd_e8, aes(y=Rate, x=Loci, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'lightgray')) +
  ylab("Percentage of True Positive") +
  guides(fill=FALSE) +
  xlab("") +
  theme_classic()
### with letters

Tposplot_e8 <- ggplot(Tpos_meansd_e8, aes(y=Rate, x=Loci, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("True Discovery Rate") +
  guides(fill="none") +
  xlab("") +
  theme(axis.text=element_text(size=12)) +geom_text(
    aes(label = c("a","c","b","a","c","b","a","b","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)



#### update the FDR :Number of false positives/Total number of all positives

FDR_e8<-data.frame(FDR=c(Fpos_pcadapt_e8/Allpositivespcadapt_e8, Fpos_rdadapt_e8/Allpositivesrda_e8,Fpos_DLmlp_e8/AllpositivesDLmlp_e8), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DNN", 100)))
FDR_mean_e8<- cbind(aggregate(FDR_e8$FDR, by=list(FDR_e8$Method), FUN = "mean"), sd=aggregate(FDR_e8$FDR, by=list(FDR_e8$Method), FUN = "sd")$x)
colnames(FDR_mean_e8)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp_e8/AllpositivesDLmlp_e8,Fpos_rdadapt_e8/Allpositivesrda_e8)
t.test(Fpos_DLmlp_e8/AllpositivesDLmlp_e8,Fpos_pcadapt_e8/Allpositivespcadapt_e8)
t.test(Fpos_rdadapt_e8/Allpositivesrda_e8,Fpos_pcadapt_e8/Allpositivespcadapt_e8)


FDRplot_e8 <- ggplot(FDR_mean_e8, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'lightgray')) +
  ylab("False Discovery Rate") +
  ylim(0,2) +  
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
### with letters
FDRplot_e8 <- ggplot(FDR_mean_e8, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,0.25) +  
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("a","b","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)

#### all 0.01
tiff(filename = "PCA_RDA_DNN_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvalue_MLPadjp-value_all_0.01_withletters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot_e8,FDRplot_e8, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()




#### FPR: Number of false positives/Total number of negatives (1000-30=970)
Fpos_e8<-data.frame(Fpos=c(Fpos_pcadapt_e8/970, Fpos_rdadapt_e8/970,Fpos_DLmlp_e8/970), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DNN", 100)))
Fpos_mean_e8<- cbind(aggregate(Fpos_e8$Fpos, by=list(Fpos_e8$Method), FUN = "mean"), sd=aggregate(Fpos_e8$Fpos, by=list(Fpos_e8$Method), FUN = "sd")$x)
colnames(Fpos_mean_e8)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp_e8/970,Fpos_rdadapt_e8/970)
t.test(Fpos_DLmlp_e8/970,Fpos_pcadapt_e8/970)
t.test(Fpos_rdadapt_e8/970,Fpos_pcadapt_e8/970)


Fposplot_e8 <- ggplot(Fpos_mean_e8, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'lightgray')) +
  ylab("False Discovery Rate") +
  ylim(0,2) +  
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
### with letters
Fposplot_e8 <- ggplot(Fpos_mean_e8, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Positive Rate") +
  ylim(0,0.025) +  
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("a","b","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)

#### all 0.01
tiff(filename = "PCA_RDA_DNN_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvalue_MLPadjp-value_all_0.01_withletters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot_e8,Fposplot_e8, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()
pdf("PCA_RDA_DNN_Scan_True_and False_discovery_rats_100simus_aPCA_RDA_qvalue_MLPdjp-value_all_0.01_withletters.pdf",width = 14,height = 7,pointsize = 12)
plot_grid(Tposplot_e8,Fposplot_e8, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()



############################ comparing arsince transformed and not PCA RDA using 0.001
#### We set p value for RDA and PCAdapt as 1e-8 to contrl FDR
Tpos_pcadapt_QTL1_e3<-unlist(lapply(p.values_pcadapt_simus, function(x) length(which(x < 1e-3 & Loci=="QT1"))))
Tpos_pcadapt_QTL2_e3<-unlist(lapply(p.values_pcadapt_simus, function(x) length(which(x < 1e-3 & Loci=="QT2"))))
Tpos_pcadapt_QTL3_e3<-unlist(lapply(p.values_pcadapt_simus, function(x) length(which(x < 1e-3 & Loci=="QT3"))))

Tpos_rdadapt_QTL1_e3<-unlist(lapply(p.values_rdadapt_simus, function(x) length(which(x < 1e-3 & Loci=="QT1"))))
Tpos_rdadapt_QTL2_e3<-unlist(lapply(p.values_rdadapt_simus, function(x) length(which(x < 1e-3 & Loci=="QT2"))))
Tpos_rdadapt_QTL3_e3<-unlist(lapply(p.values_rdadapt_simus, function(x) length(which(x < 1e-3 & Loci=="QT3"))))

Tpos_DLmlp_QTL1_e10<-unlist(lapply(p.values_DLmlp_simusarcsin, function(x) length(which(x < 1e-10  & Loci=="QT1"))))
Tpos_DLmlp_QTL2_e10<-unlist(lapply(p.values_DLmlp_simusarcsin, function(x) length(which(x < 1e-10 & Loci=="QT2"))))
Tpos_DLmlp_QTL3_e10<-unlist(lapply(p.values_DLmlp_simusarcsin, function(x) length(which(x < 1e-10  & Loci=="QT3"))))

## adju pvalue, do these later to compare the p.values and then choose one form


## False postives 
Fpos_pcadapt_e3<-unlist(lapply(p.values_pcadapt_simus, function(x) length(which(x < 1e-3 & Loci=="Neutral"))))
Fpos_rdadapt_e3<-unlist(lapply(p.values_rdadapt_simus, function(x) length(which(x < 1e-3 & Loci=="Neutral"))))

Fpos_DLmlp_e10<-unlist(lapply(p.values_DLmlp_simusarcsin, function(x) length(which(x <1e-10 & Loci=="Neutral"))))


### t test to see the difference
t.test(Tpos_DLmlp_QTL1_e10/10,Tpos_rdadapt_QTL1_e3/10)
t.test(Tpos_DLmlp_QTL1_e10/10,Tpos_pcadapt_QTL1_e3/10)
t.test(Tpos_rdadapt_QTL1_e3/10,Tpos_pcadapt_QTL1_e3/10)


t.test(Tpos_DLmlp_QTL2_e10/10,Tpos_rdadapt_QTL2_e3/10)
t.test(Tpos_DLmlp_QTL2_e10/10,Tpos_pcadapt_QTL2_e3/10)
t.test(Tpos_rdadapt_QTL2_e3/10,Tpos_pcadapt_QTL2_e3/10)

t.test(Tpos_DLmlp_QTL3_e10/10,Tpos_rdadapt_QTL3_e3/10)
t.test(Tpos_DLmlp_QTL3_e10/10,Tpos_pcadapt_QTL3_e3/10)
t.test(Tpos_rdadapt_QTL3_e3/10,Tpos_pcadapt_QTL3_e3/10)


Allpositivespcadapt_e3=unlist(lapply(p.values_pcadapt_simus, function(x) length(which(x < 1e-3))))
Allpositivesrda_e3=unlist(lapply(p.values_rdadapt_simus, function(x) length(which(x < 1e-3))))
AllpositivesDLmlp_e10=unlist(lapply(p.values_DLmlp_simusarcsin, function(x) length(which(x < 1e-10))))


### True Discovery rate
Tpos_e10_3<-data.frame(Tpos=c(Tpos_pcadapt_QTL1_e3/10, Tpos_pcadapt_QTL2_e3/10, Tpos_pcadapt_QTL3_e3/10, Tpos_rdadapt_QTL1_e3/10, Tpos_rdadapt_QTL2_e3/10, Tpos_rdadapt_QTL3_e3/10,Tpos_DLmlp_QTL1_e10/10, Tpos_DLmlp_QTL2_e10/10, Tpos_DLmlp_QTL3_e10/10), Method=c(rep("PCA", 300), rep("RDA", 300),rep("DNN", 300)), Loci=c(rep("QT1", 100), rep("QT2", 100), rep("QT3", 100), rep("QT1", 100), rep("QT2", 100), rep("QT3", 100),rep("QT1", 100), rep("QT2", 100), rep("QT3", 100)))
Tpos_meansd_e10_3<-cbind(aggregate(Tpos_e10_3$Tpos, by=list(Tpos_e10_3$Method,Tpos_e10_3$Loci), FUN = "mean"), sd=aggregate(Tpos_e10_3$Tpos, by=list(Tpos_e10_3$Method,Tpos_e10_3$Loci), FUN = "sd")$x)
colnames(Tpos_meansd_e10_3)[1:3]<-c("Method", "Loci", "Rate")




## Visualisation

### with letters

Tposplot_e10_3 <- ggplot(Tpos_meansd_e10_3, aes(y=Rate, x=Loci, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("True Discovery Rate") +
  guides(fill="none") +
  xlab("") +
  theme(axis.text=element_text(size=12)) +geom_text(
    aes(label = c("a","c","b","a","b","a","a","c","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)



#### update the FDR :Number of false positives/Total number of all positives

FDR_e10_3<-data.frame(FDR=c(Fpos_pcadapt_e3/Allpositivespcadapt_e3, Fpos_rdadapt_e3/Allpositivesrda_e3,Fpos_DLmlp_e10/AllpositivesDLmlp_e10), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DNN", 100)))
FDR_mean_e10_3<- cbind(aggregate(FDR_e10_3$FDR, by=list(FDR_e10_3$Method), FUN = "mean"), sd=aggregate(FDR_e10_3$FDR, by=list(FDR_e10_3$Method), FUN = "sd")$x)
colnames(FDR_mean_e10_3)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp_e10/AllpositivesDLmlp_e10,Fpos_rdadapt_e3/Allpositivesrda_e3)
t.test(Fpos_DLmlp_e10/AllpositivesDLmlp_e10,Fpos_pcadapt_e3/Allpositivespcadapt_e3)
t.test(Fpos_rdadapt_e3/Allpositivesrda_e3,Fpos_pcadapt_e3/Allpositivespcadapt_e3)


FDRplot_e10_3 <- ggplot(FDR_mean_e10_3, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'lightgray')) +
  ylab("False Discovery Rate") +
  ylim(0,2) +  
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

### with letters
FDRplot_e10_3 <- ggplot(FDR_mean_e10_3, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,0.5) +  
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("b","a","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)

#### all 0.01
tiff(filename = "PCA_RDA_DNN_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvalue_MLPadjp-value_all_0.01_withletters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot_e10_3,FDRplot_e10_3, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()




#### FPR: Number of false positives/Total number of negatives (1000-30=970)
Fpos_e10_3<-data.frame(Fpos=c(Fpos_pcadapt_e3/970, Fpos_rdadapt_e3/970,Fpos_DLmlp_e10/970), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DNN", 100)))
Fpos_mean_e10_3<- cbind(aggregate(Fpos_e10_3$Fpos, by=list(Fpos_e10_3$Method), FUN = "mean"), sd=aggregate(Fpos_e10_3$Fpos, by=list(Fpos_e10_3$Method), FUN = "sd")$x)
colnames(Fpos_mean_e10_3)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp_e10/970,Fpos_rdadapt_e3/970) ### BF correction/1000
t.test(Fpos_DLmlp_e10/970,Fpos_pcadapt_e3/970) 
t.test(Fpos_rdadapt_e3/970,Fpos_pcadapt_e3/970)


Fposplot_e10_3<- ggplot(Fpos_mean_e10_3, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'lightgray')) +
  ylab("False Discovery Rate") +
  ylim(0,2) +  
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
### with letters
Fposplot_e10_3 <- ggplot(Fpos_mean_e10_3, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Positive Rate") +
  ylim(0,0.025) +  
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("b","a","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)

#### all 0.01
tiff(filename = "PCA_RDA_DNN_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvalue_MLPadjp-value_all_0.01_withletters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot_e10_3,Fposplot_e10_3, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()
pdf("PCA_RDA_DNN_Scan_True_and False_discovery_rats_100simus_aPCA_RDA_qvalue_MLPdjp-value_all_0.01_withletters.pdf",width = 14,height = 7,pointsize = 12)
plot_grid(Tposplot_e10_3,Fposplot_e10_3, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()

############################################################################################################################################
################ e-8: DL; e-3:PCA,RDA
### t test to see the difference
t.test(Tpos_DLmlp_QTL1_e8/10,Tpos_rdadapt_QTL1_e3/10)
t.test(Tpos_DLmlp_QTL1_e8/10,Tpos_pcadapt_QTL1_e3/10)
t.test(Tpos_rdadapt_QTL1_e3/10,Tpos_pcadapt_QTL1_e3/10)


t.test(Tpos_DLmlp_QTL2_e8/10,Tpos_rdadapt_QTL2_e3/10)
t.test(Tpos_DLmlp_QTL2_e8/10,Tpos_pcadapt_QTL2_e3/10)
t.test(Tpos_rdadapt_QTL2_e3/10,Tpos_pcadapt_QTL2_e3/10)

t.test(Tpos_DLmlp_QTL3_e8/10,Tpos_rdadapt_QTL3_e3/10)
t.test(Tpos_DLmlp_QTL3_e8/10,Tpos_pcadapt_QTL3_e3/10)
t.test(Tpos_rdadapt_QTL3_e3/10,Tpos_pcadapt_QTL3_e3/10)


Allpositivespcadapt_e3=unlist(lapply(p.values_pcadapt_simus, function(x) length(which(x < 1e-3))))
Allpositivesrda_e3=unlist(lapply(p.values_rdadapt_simus, function(x) length(which(x < 1e-3))))
AllpositivesDLmlp_e8=unlist(lapply(p.values_DLmlp_simusarcsin, function(x) length(which(x < 1e-8))))


### True Discovery rate
Tpos_e8_3<-data.frame(Tpos=c(Tpos_pcadapt_QTL1_e3/10, Tpos_pcadapt_QTL2_e3/10, Tpos_pcadapt_QTL3_e3/10, Tpos_rdadapt_QTL1_e3/10, Tpos_rdadapt_QTL2_e3/10, Tpos_rdadapt_QTL3_e3/10,Tpos_DLmlp_QTL1_e8/10, Tpos_DLmlp_QTL2_e8/10, Tpos_DLmlp_QTL3_e8/10), Method=c(rep("PCA", 300), rep("RDA", 300),rep("DNN", 300)), Loci=c(rep("QT1", 100), rep("QT2", 100), rep("QT3", 100), rep("QT1", 100), rep("QT2", 100), rep("QT3", 100),rep("QT1", 100), rep("QT2", 100), rep("QT3", 100)))
Tpos_meansd_e8_3<-cbind(aggregate(Tpos_e8_3$Tpos, by=list(Tpos_e8_3$Method,Tpos_e8_3$Loci), FUN = "mean"), sd=aggregate(Tpos_e8_3$Tpos, by=list(Tpos_e8_3$Method,Tpos_e8_3$Loci), FUN = "sd")$x)
colnames(Tpos_meansd_e8_3)[1:3]<-c("Method", "Loci", "Rate")




## Visualisation

### with letters

Tposplot_e8_3 <- ggplot(Tpos_meansd_e8_3, aes(y=Rate, x=Loci, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("True Discovery Rate") +
  guides(fill="none") +
  xlab("") +
  theme(axis.text=element_text(size=12)) +geom_text(
    aes(label = c("a","c","b","a","c","b","a","c","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)



#### update the FDR :Number of false positives/Total number of all positives

FDR_e8_3<-data.frame(FDR=c(Fpos_pcadapt_e3/Allpositivespcadapt_e3, Fpos_rdadapt_e3/Allpositivesrda_e3,Fpos_DLmlp_e8/AllpositivesDLmlp_e8), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DNN", 100)))
FDR_mean_e8_3<- cbind(aggregate(FDR_e8_3$FDR, by=list(FDR_e8_3$Method), FUN = "mean"), sd=aggregate(FDR_e8_3$FDR, by=list(FDR_e8_3$Method), FUN = "sd")$x)
colnames(FDR_mean_e8_3)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp_e8/AllpositivesDLmlp_e8,Fpos_rdadapt_e3/Allpositivesrda_e3)
t.test(Fpos_DLmlp_e8/AllpositivesDLmlp_e8,Fpos_pcadapt_e3/Allpositivespcadapt_e3)
t.test(Fpos_rdadapt_e3/Allpositivesrda_e3,Fpos_pcadapt_e3/Allpositivespcadapt_e3)


FDRplot_e8_3 <- ggplot(FDR_mean_e8_3, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'lightgray')) +
  ylab("False Discovery Rate") +
  ylim(0,2) +  
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

### with letters
FDRplot_e8_3 <- ggplot(FDR_mean_e8_3, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,0.5) +  
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("b","a","c"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)

#### all 0.01
tiff(filename = "PCA_RDA_DNN_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvalue_MLPadjp-value_all_0.01_withletters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot_e8_3,FDRplot_e8_3, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()




#### FPR: Number of false positives/Total number of negatives (1000-30=970)
Fpos_e8_3<-data.frame(Fpos=c(Fpos_pcadapt_e3/970, Fpos_rdadapt_e3/970,Fpos_DLmlp_e8/970), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DNN", 100)))
Fpos_mean_e8_3<- cbind(aggregate(Fpos_e8_3$Fpos, by=list(Fpos_e8_3$Method), FUN = "mean"), sd=aggregate(Fpos_e8_3$Fpos, by=list(Fpos_e8_3$Method), FUN = "sd")$x)
colnames(Fpos_mean_e8_3)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp_e8/970,Fpos_rdadapt_e3/970)
t.test(Fpos_DLmlp_e8/970,Fpos_pcadapt_e3/970)
t.test(Fpos_rdadapt_e3/970,Fpos_pcadapt_e3/970)


Fposplot_e8_3<- ggplot(Fpos_mean_e8_3, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'lightgray')) +
  ylab("False Discovery Rate") +
  ylim(0,2) +  
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
### with letters
Fposplot_e8_3 <- ggplot(Fpos_mean_e8_3, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Positive Rate") +
  ylim(0,0.025) +  
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("b","a","c"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)

#### all 0.01
tiff(filename = "PCA_RDA_DNN_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvalue_MLPadjp-value_all_0.01_withletters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot_e8_3,Fposplot_e8_3, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()
pdf("PCA_RDA_DNN_Scan_True_and False_discovery_rats_100simus_aPCA_RDA_qvalue_MLPdjp-value_all_0.01_withletters.pdf",width = 14,height = 7,pointsize = 12)
plot_grid(Tposplot_e8_3,Fposplot_e8_3, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()
###########20-04-2022


##########################  e-12

Tpos_DLmlp_QTL1_e12<-unlist(lapply(p.values_DLmlp_simusarcsin, function(x) length(which(x < 1e-12  & Loci=="QT1"))))
Tpos_DLmlp_QTL2_e12<-unlist(lapply(p.values_DLmlp_simusarcsin, function(x) length(which(x < 1e-12 & Loci=="QT2"))))
Tpos_DLmlp_QTL3_e12<-unlist(lapply(p.values_DLmlp_simusarcsin, function(x) length(which(x < 1e-12  & Loci=="QT3"))))

## adju pvalue, do these later to compare the p.values and then choose one form


## False postives 
Fpos_pcadapt_e3<-unlist(lapply(p.values_pcadapt_simus, function(x) length(which(x < 1e-3 & Loci=="Neutral"))))
Fpos_rdadapt_e3<-unlist(lapply(p.values_rdadapt_simus, function(x) length(which(x < 1e-3 & Loci=="Neutral"))))

Fpos_DLmlp_e12<-unlist(lapply(p.values_DLmlp_simusarcsin, function(x) length(which(x <1e-12 & Loci=="Neutral"))))


### t test to see the difference
t.test(Tpos_DLmlp_QTL1_e12/10,Tpos_rdadapt_QTL1_e3/10)
t.test(Tpos_DLmlp_QTL1_e12/10,Tpos_pcadapt_QTL1_e3/10)
t.test(Tpos_rdadapt_QTL1_e3/10,Tpos_pcadapt_QTL1_e3/10)


t.test(Tpos_DLmlp_QTL2_e12/10,Tpos_rdadapt_QTL2_e3/10)
t.test(Tpos_DLmlp_QTL2_e12/10,Tpos_pcadapt_QTL2_e3/10)
t.test(Tpos_rdadapt_QTL2_e3/10,Tpos_pcadapt_QTL2_e3/10)

t.test(Tpos_DLmlp_QTL3_e12/10,Tpos_rdadapt_QTL3_e3/10)
t.test(Tpos_DLmlp_QTL3_e12/10,Tpos_pcadapt_QTL3_e3/10)
t.test(Tpos_rdadapt_QTL3_e3/10,Tpos_pcadapt_QTL3_e3/10)


Allpositivespcadapt_e3=unlist(lapply(p.values_pcadapt_simus, function(x) length(which(x < 1e-3))))
Allpositivesrda_e3=unlist(lapply(p.values_rdadapt_simus, function(x) length(which(x < 1e-3))))
AllpositivesDLmlp_e12=unlist(lapply(p.values_DLmlp_simusarcsin, function(x) length(which(x < 1e-12))))


### True Discovery rate
Tpos_e12_3<-data.frame(Tpos=c(Tpos_pcadapt_QTL1_e3/10, Tpos_pcadapt_QTL2_e3/10, Tpos_pcadapt_QTL3_e3/10, Tpos_rdadapt_QTL1_e3/10, Tpos_rdadapt_QTL2_e3/10, Tpos_rdadapt_QTL3_e3/10,Tpos_DLmlp_QTL1_e12/10, Tpos_DLmlp_QTL2_e12/10, Tpos_DLmlp_QTL3_e12/10), Method=c(rep("PCA", 300), rep("RDA", 300),rep("DNN", 300)), Loci=c(rep("QT1", 100), rep("QT2", 100), rep("QT3", 100), rep("QT1", 100), rep("QT2", 100), rep("QT3", 100),rep("QT1", 100), rep("QT2", 100), rep("QT3", 100)))
Tpos_meansd_e12_3<-cbind(aggregate(Tpos_e12_3$Tpos, by=list(Tpos_e12_3$Method,Tpos_e12_3$Loci), FUN = "mean"), sd=aggregate(Tpos_e12_3$Tpos, by=list(Tpos_e12_3$Method,Tpos_e12_3$Loci), FUN = "sd")$x)
colnames(Tpos_meansd_e12_3)[1:3]<-c("Method", "Loci", "Rate")




## Visualisation

### with letters

Tposplot_e12_3 <- ggplot(Tpos_meansd_e12_3, aes(y=Rate, x=Loci, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("True Discovery Rate") +
  guides(fill="none") +
  xlab("") +
  theme(axis.text=element_text(size=12)) +geom_text(
    aes(label = c("a","c","b","a","b","a","a","c","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)



#### update the FDR :Number of false positives/Total number of all positives

FDR_e12_3<-data.frame(FDR=c(Fpos_pcadapt_e3/Allpositivespcadapt_e3, Fpos_rdadapt_e3/Allpositivesrda_e3,Fpos_DLmlp_e12/AllpositivesDLmlp_e12), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DNN", 100)))
FDR_mean_e12_3<- cbind(aggregate(FDR_e12_3$FDR, by=list(FDR_e12_3$Method), FUN = "mean"), sd=aggregate(FDR_e12_3$FDR, by=list(FDR_e12_3$Method), FUN = "sd")$x)
colnames(FDR_mean_e12_3)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp_e12/AllpositivesDLmlp_e12,Fpos_rdadapt_e3/Allpositivesrda_e3)
t.test(Fpos_DLmlp_e12/AllpositivesDLmlp_e12,Fpos_pcadapt_e3/Allpositivespcadapt_e3)
t.test(Fpos_rdadapt_e3/Allpositivesrda_e3,Fpos_pcadapt_e3/Allpositivespcadapt_e3)


FDRplot_e12_3 <- ggplot(FDR_mean_e12_3, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'lightgray')) +
  ylab("False Discovery Rate") +
  ylim(0,2) +  
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

### with letters
FDRplot_e12_3 <- ggplot(FDR_mean_e12_3, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,0.5) +  
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("b","a","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)

#### all 0.01
tiff(filename = "PCA_RDA_DNN_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvalue_MLPadjp-value_all_0.01_withletters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot_e12_3,FDRplot_e12_3, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()




#### FPR: Number of false positives/Total number of negatives (1000-30=970)
Fpos_e12_3<-data.frame(Fpos=c(Fpos_pcadapt_e3/970, Fpos_rdadapt_e3/970,Fpos_DLmlp_e12/970), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DNN", 100)))
Fpos_mean_e12_3<- cbind(aggregate(Fpos_e12_3$Fpos, by=list(Fpos_e12_3$Method), FUN = "mean"), sd=aggregate(Fpos_e12_3$Fpos, by=list(Fpos_e12_3$Method), FUN = "sd")$x)
colnames(Fpos_mean_e12_3)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp_e12/970,Fpos_rdadapt_e3/970) ### BF correction/1000
t.test(Fpos_DLmlp_e12/970,Fpos_pcadapt_e3/970) 
t.test(Fpos_rdadapt_e3/970,Fpos_pcadapt_e3/970)


Fposplot_e12_3<- ggplot(Fpos_mean_e12_3, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'lightgray')) +
  ylab("False Discovery Rate") +
  ylim(0,2) +  
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
### with letters
Fposplot_e12_3 <- ggplot(Fpos_mean_e12_3, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Positive Rate") +
  ylim(0,0.025) +  
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("b","a","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)

#### all 0.01
tiff(filename = "PCA_RDA_DNN_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvalue_MLPadjp-value_all_0.01_withletters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot_e12_3,Fposplot_e12_3, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()
pdf("PCA_RDA_DNN_Scan_True_and False_discovery_rats_100simus_aPCA_RDA_qvalue_MLPdjp-value_all_0.01_withletters.pdf",width = 14,height = 7,pointsize = 12)
plot_grid(Tposplot_e12_3,Fposplot_e12_3, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()


#####################################################combing plots###################################################################
Tposplot_e8 <- ggplot(Tpos_meansd_e8, aes(y=Rate, x=Loci, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("True Discovery Rate") +
  guides(fill="none") +
  xlab("") +
  theme(axis.text=element_text(size=12)) +geom_text(
    aes(label = c("a","c","b","a","c","b","a","b","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)

FDRplot_e8 <- ggplot(FDR_mean_e8, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,0.5) +  
  guides(fill="none") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("a","b","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)

Fposplot_e8 <- ggplot(Fpos_mean_e8, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Positive Rate") +
  ylim(0,0.01) +  
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("a","b","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)


############################e-8, e-3

Tposplot_e8_3 <- ggplot(Tpos_meansd_e8_3, aes(y=Rate, x=Loci, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("True Discovery Rate") +
  guides(fill="none") +
  xlab("") +
  theme(axis.text=element_text(size=12)) +geom_text(
    aes(label = c("a","c","b","a","c","b","a","c","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)

FDRplot_e8_3 <- ggplot(FDR_mean_e8_3, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,0.5) +  
  guides(fill="none") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("b","a","c"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)

Fposplot_e8_3 <- ggplot(Fpos_mean_e8_3, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Positive Rate") +
  ylim(0,0.01) +  
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("b","a","c"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)


############################ e-10, e-3

Tposplot_e10_3 <- ggplot(Tpos_meansd_e10_3, aes(y=Rate, x=Loci, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("True Discovery Rate") +
  guides(fill="none") +
  xlab("") +
  theme(axis.text=element_text(size=12)) +geom_text(
    aes(label = c("a","c","b","a","b","a","a","c","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)

FDRplot_e10_3 <- ggplot(FDR_mean_e10_3, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,0.5) +  
  guides(fill="none") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("b","a","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)

Fposplot_e10_3 <- ggplot(Fpos_mean_e10_3, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Positive Rate") +
  ylim(0,0.01) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("b","a","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)
plot_grid(Tposplot_e8,FDRplot_e8,Fposplot_e8,Tposplot_e8_3,FDRplot_e8_3,Fposplot_e8_3,Tposplot_e10_3,FDRplot_e10_3,Fposplot_e10_3, labels = "AUTO", align = "hv", rel_widths = c(1, .7))

#################################






##### combining figures


tiff(filename = "PCA_RDA_DL_MLP_Scan_power_1_withletters20220420.tif",width = 18,height = 18,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot_e8,FDRplot_e8,Fposplot_e8,Tposplot_e8_3,FDRplot_e8_3,Fposplot_e8_3,Tposplot_e10_3,FDRplot_e10_3,Fposplot_e10_3, ncol =3, nrow=3,labels = c("A","B","C","D","E","F","G","H"),font.label = list(size = 20),align = "h" )
dev.off()

pdf("PCA_RDA_DL_MLP_Scan_power_1_withletters20220420.pdf",width = 18,height = 18,pointsize = 12)
plot_grid(Tposplot_e8,FDRplot_e8,Fposplot_e8,Tposplot_e8_3,FDRplot_e8_3,Fposplot_e8_3,Tposplot_e10_3,FDRplot_e10_3,Fposplot_e10_3,ncol =3, nrow=3,labels =c("A","B","C","D","E","F","G","H"), font.label = list(size = 20),align = "h")
dev.off()
#### small letters
tiff(filename = "PCA_RDA_DL_MLP_Scan_power_1_withletters_small20220420.tif",width = 18,height = 18,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot_e8,FDRplot_e8,Fposplot_e8,Tposplot_e8_3,FDRplot_e8_3,Fposplot_e8_3,Tposplot_e10_3,FDRplot_e10_3,Fposplot_e10_3,ncol =3, nrow=3,labels = c("a","b","c","d","e","f","g","h"),font.label = list(size = 20), align = "h")
dev.off()

pdf("PCA_RDA_DL_MLP_Scan_power_1_withletters_small20220420.pdf",width = 18,height = 18,pointsize = 12)
plot_grid(Tposplot_e8,FDRplot_e8,Fposplot_e8,Tposplot_e8_3,FDRplot_e8_3,Fposplot_e8_3,Tposplot_e10_3,FDRplot_e10_3,Fposplot_e10_3,ncol =3, nrow=3,labels = c("a","b","c","d","e","f","g","h"),font.label = list(size = 20), align = "h")
dev.off()

##############
tiff(filename = "PCA_RDA_DL_MLP_Scan_power__withletters20220420.tif",width = 18,height = 18,units = "in",pointsize = 12,res=600,compression = "zip+p")
ggarrange(Tposplot_e8,FDRplot_e8,Fposplot_e8,Tposplot_e8_3,FDRplot_e8_3,Fposplot_e8_3,Tposplot_e10_3,FDRplot_e10_3,Fposplot_e10_3, ncol =3, nrow=3,labels = c("A","B","C","D","E","F","G","H"),font.label = list(size = 20),align = "h" )
dev.off()

pdf("PCA_RDA_DL_MLP_Scan_power__withletters20220420.pdf",width = 18,height = 18,pointsize = 12)
ggarrange(Tposplot_e8,FDRplot_e8,Fposplot_e8,Tposplot_e8_3,FDRplot_e8_3,Fposplot_e8_3,Tposplot_e10_3,FDRplot_e10_3,Fposplot_e10_3,ncol =3, nrow=3,labels =c("A","B","C","D","E","F","G","H"), font.label = list(size = 20),align = "h")
dev.off()
#### small letters
tiff(filename = "PCA_RDA_DL_MLP_Scan_power__withletters_small20220420.tif",width = 18,height = 18,units = "in",pointsize = 12,res=600,compression = "zip+p")
ggarrange(Tposplot_e8,FDRplot_e8,Fposplot_e8,Tposplot_e8_3,FDRplot_e8_3,Fposplot_e8_3,Tposplot_e10_3,FDRplot_e10_3,Fposplot_e10_3,ncol =3, nrow=3,labels = c("a","b","c","d","e","f","g","h"), font.label = list(size = 20),align = "h")
dev.off()

pdf("PCA_RDA_DL_MLP_Scan_power__withletters_small20220420.pdf",width = 18,height = 18,pointsize = 12)
ggarrange(Tposplot_e8,FDRplot_e8,Fposplot_e8,Tposplot_e8_3,FDRplot_e8_3,Fposplot_e8_3,Tposplot_e10_3,FDRplot_e10_3,Fposplot_e10_3,ncol =3, nrow=3,labels = c("a","b","c","d","e","f","g","h"),font.label = list(size = 20), align = "h")
dev.off()







