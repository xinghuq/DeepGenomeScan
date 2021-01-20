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
impfile <- list.files(impdir, pattern="*_mlpneuralnet_importance.csv", full.names = TRUE)
## 100 simulations and 10 envir variables, 

impSim_all=lapply(impfile, read.csv) ### including all 100 simulations and 10 variables each simuilations, 
impfile[11:20]## check, each 10 file is a simulation
### x is the loci ID and Overall is the overall importance, 1000 list: including 100 simulations and 10 envir

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

DLsim1=DLdata[111:120] ## This is the simulation 1, check
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
chunk(1:1000,100)
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
chunk2(1:1000,100)
samplesim=chunk2(1:1000,100)


DLqvalues<-function(DL_data,K)
{
  loadings<-DL_data# [,1:as.numeric(K)]
  resscale <- apply(loadings, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="donostah")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_DL<-qval$qvalues
  padj <- p.adjust(reschi2test,method="bonferroni")
  return(data.frame(p.values=reschi2test, q.values=q.values_DL,padj=padj))
}
###The covRob function selects a robust covariance estimator that is likely to provide a good estimate in a reasonable amount of time. Presently this selection is based on the problem size. The Donoho-Stahel estimator is used if there are less than 1000 observations and less than 10 variables or less than 5000 observations and less than 5 variables. If there are less than 50000 observations and less than 20 variables then the MCD is used. For larger problems, the Orthogonalized Quadrant Correlation estimator is used.
Simqvalue=DLqvalues(DLsim1,10)

## Manhattan plot
ggplot() +
  geom_point(aes(x=which(Loci=="Neutral"), y=-log10(Simqvalue[-which(Loci!="Neutral"),2])), col = "gray83") +
  geom_point(aes(x=which(Loci!="Neutral"), y=-log10(Simqvalue[-which(Loci=="Neutral"),2]), colour = Selected_Loci)) +
  xlab("SNPs") + ylab("DL with condition)") +ylim(c(0,300))+theme_bw()



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

DLqvalues<-function(DL_data,K)
{
  DL_importance<-DL_data  #[,1:as.numeric(nenvir)]
  resscale <- apply(DL_importance, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="donostah")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_DL<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_DL))
}



#donostah is good if variables is less than 1000

Loci<-rep("Neutral", 1000)
Loci[c(1,11,21,31,41,51,61,71,81,91)]<-"QTL1"
Loci[c(101,111,121,131,141,151,161,171,181,191)]<-"QTL2"
Loci[c(201,211,221,231,241,251,261,271,281,291)]<-"QTL3"
### q values
p.values<-data.frame(index = c(1:1000), pcadapt = -log10(qvalue(pcadapt$pvalues)$qvalues), rdadapt = -log10(res_rdadapt[,2]), DL_mlp=-log10(Simqvalue$q.values),Loci=Loci)
p.values<-data.frame(index = c(1:1000), pcadapt = -log10(qvalue(pcadapt$pvalues)$qvalues), rdadapt = -log10(res_rdadapt[,2]), DL_mlp=-log10(Simqvalue$padj),Loci=Loci)

p.values<-data.frame(index = c(1:1000), pcadapt = -log10(pcadapt$pvalues), rdadapt = -log10(res_rdadapt[,1]), DL_mlp=-log10(Simqvalue$q.values),Loci=Loci)

Selected_Loci<-p.values$Loci[-which(p.values$Loci=="Neutral")]

#### Fig. S3


S3p1 <- ggplot() +
  geom_point(aes(x=p.values$index[-which(p.values$Loci!="Neutral")], y=p.values$pcadapt[-which(p.values$Loci!="Neutral")]), col = "gray83") +
  geom_point(aes(x=as.vector(p.values$index[-which(p.values$Loci=="Neutral")]), y=as.vector(p.values$pcadapt[-which(p.values$Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNP (with maf > 0.05)") + ylab("PCADAPT -log10(q-values)") +
  ylim(0,6) +
  theme_bw()
S3p2 <- ggplot() +
  geom_point(aes(x=p.values$index[-which(p.values$Loci!="Neutral")], y=p.values$rdadapt[-which(p.values$Loci!="Neutral")]), col = "gray83") +
  geom_point(aes(x=as.vector(p.values$index[-which(p.values$Loci=="Neutral")]), y=as.vector(p.values$rdadapt[-which(p.values$Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNP (with maf > 0.05)") + ylab("RDA_Scan -log10(q-values)") +
  ylim(0,6) +
  theme_bw()
S3p3 <- ggplot() +
  geom_point(aes(x=p.values$index[-which(p.values$Loci!="Neutral")], y=p.values$DL_mlp[-which(p.values$Loci!="Neutral")]), col = "gray83") +
  geom_point(aes(x=as.vector(p.values$index[-which(p.values$Loci=="Neutral")]), y=as.vector(p.values$DL_mlp[-which(p.values$Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNP (with maf > 0.05)") + ylab("DL_MLP -log10(q-values)") +
  ylim(0,300) +
  theme_bw()

ggarrange(S3p1, S3p2, S3p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE, legend="right")

tiff(filename = "PCA_RDA_DL_MLP_Scan_Sim1_qvalue_all1_small_letters.tif",width = 21,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
ggarrange(S3p1, S3p2, S3p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE, font.label = list(size=18),  legend="right")
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_Sim1_qvalue_all1_small_letters.pdf",width = 21,height = 7,pointsize = 12)
ggarrange(S3p1, S3p2, S3p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE, font.label = list(size=18),  legend="right")
dev.off()

tiff(filename = "PCA_RDA_DL_MLP_Scan_Sim1_14_7_qvalue_all1_small_letters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
ggarrange(S3p1, S3p2, S3p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE,font.label = list(size=18),   legend="right")
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_Sim1_14_7_qvalue_all1_small_letters.pdf",width = 14,height = 7,pointsize = 12)
ggarrange(S3p1, S3p2, S3p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE,font.label = list(size=18),   legend="right")
dev.off()



p1 <- ggplot() +
  geom_point(aes(x=p.values$index[-which(p.values$Loci!="Neutral")], y=p.values$pcadapt[-which(p.values$Loci!="Neutral")]), col = "gray83") +
  geom_point(aes(x=as.vector(p.values$index[-which(p.values$Loci=="Neutral")]), y=as.vector(p.values$pcadapt[-which(p.values$Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNP (with maf > 0.05)") + ylab("PCADAPT -log10(q-values)") +
  ylim(0,300) +
  theme_bw()
p2 <- ggplot() +
  geom_point(aes(x=p.values$index[-which(p.values$Loci!="Neutral")], y=p.values$rdadapt[-which(p.values$Loci!="Neutral")]), col = "gray83") +
  geom_point(aes(x=as.vector(p.values$index[-which(p.values$Loci=="Neutral")]), y=as.vector(p.values$rdadapt[-which(p.values$Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNP (with maf > 0.05)") + ylab("RDA_Scan -log10(q-values)") +
  ylim(0,300) +
  theme_bw()
p3 <- ggplot() +
  geom_point(aes(x=p.values$index[-which(p.values$Loci!="Neutral")], y=p.values$DL_mlp[-which(p.values$Loci!="Neutral")]), col = "gray83") +
  geom_point(aes(x=as.vector(p.values$index[-which(p.values$Loci=="Neutral")]), y=as.vector(p.values$DL_mlp[-which(p.values$Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNP (with maf > 0.05)") + ylab("DL_MLP -log10(q-values)") +
  ylim(0,300) +
  theme_bw()

#### 05-09-2020 update change front letters as small letters


tiff(filename = "PCA_RDA_DL_MLP_Scan_Sim1_log_qvalues1_small_letters.tif",width = 21,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
ggarrange(p1, p2, p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE,font.label = list(size=18), legend="right")
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_Sim1_log_qvalues1_small_letters.pdf",width = 21,height = 7,pointsize = 12)
ggarrange(p1, p2, p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE, font.label = list(size=18), legend="right")
dev.off()

tiff(filename = "PCA_RDA_DL_MLP_Scan_Sim1_14_7_log_qvalues1_small_letters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
ggarrange(p1, p2, p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE,font.label = list(size=18),  legend="right")
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_Sim1_14_7_log_qvalues1_small_letters.pdf",width = 14,height = 7,pointsize = 12)
ggarrange(p1, p2, p3,labels = c("a", "b","c"), ncol=3, common.legend = TRUE, font.label = list(size=18), legend="right")
dev.off()
########
save.image(file = "SimData_05_09_2020.RData")
##27-10-2020 update, produce additional figures re-aggrange the scale to macth q-value=50 as a new balance
#### Fig. S3


S3p11 <- ggplot() +
  geom_point(aes(x=p.values$index[-which(p.values$Loci!="Neutral")], y=p.values$pcadapt[-which(p.values$Loci!="Neutral")]), col = "gray83") +
  geom_point(aes(x=as.vector(p.values$index[-which(p.values$Loci=="Neutral")]), y=as.vector(p.values$pcadapt[-which(p.values$Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNP (with maf > 0.05)") + ylab("PCADAPT -log10(q-values)") +
  ylim(0,8) +
  theme_bw()
S3p21 <- ggplot() +
  geom_point(aes(x=p.values$index[-which(p.values$Loci!="Neutral")], y=p.values$rdadapt[-which(p.values$Loci!="Neutral")]), col = "gray83") +
  geom_point(aes(x=as.vector(p.values$index[-which(p.values$Loci=="Neutral")]), y=as.vector(p.values$rdadapt[-which(p.values$Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNP (with maf > 0.05)") + ylab("RDA_Scan -log10(q-values)") +
  ylim(0,8) +
  theme_bw()
S3p31 <- ggplot() +
  geom_point(aes(x=p.values$index[-which(p.values$Loci!="Neutral")], y=p.values$DL_mlp[-which(p.values$Loci!="Neutral")]), col = "gray83") +
  geom_point(aes(x=as.vector(p.values$index[-which(p.values$Loci=="Neutral")]), y=as.vector(p.values$DL_mlp[-which(p.values$Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNP (with maf > 0.05)") + ylab("DL_MLP -log10(q-values)") +
  ylim(0,200) +
  theme_bw()


ggarrange(S3p11, S3p21, S3p31,labels = c("a", "b","c"), ncol=3, common.legend = TRUE, legend="right")

tiff(filename = "PCA_RDA_DL_MLP_Scan_Sim1_qvalue_0.01vs50_small_letters.tif",width = 21,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
ggarrange(S3p11, S3p21, S3p31,labels = c("a", "b","c"), ncol=3, common.legend = TRUE, font.label = list(size=18),  legend="right")
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_Sim1_qvalue_0.01vs50_small_letters.pdf",width = 21,height = 7,pointsize = 12)
ggarrange(S3p11, S3p21, S3p31,labels = c("a", "b","c"), ncol=3, common.legend = TRUE, font.label = list(size=18),  legend="right")
dev.off()

tiff(filename = "PCA_RDA_DL_MLP_Scan_Sim1_14_7_qvalue_0.01vs50_small_letters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
ggarrange(S3p11, S3p21, S3p31,labels = c("a", "b","c"), ncol=3, common.legend = TRUE,font.label = list(size=18),   legend="right")
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_Sim1_14_7_qvalue_0.01vs50_small_letters.pdf",width = 14,height = 7,pointsize = 12)
ggarrange(S3p11, S3p21, S3p31,labels = c("a", "b","c"), ncol=3, common.legend = TRUE,font.label = list(size=18),   legend="right")
dev.off()

#### upletters

tiff(filename = "PCA_RDA_DL_MLP_Scan_Sim1_qvalue_0.01vs50_initial_letters.tif",width = 21,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
ggarrange(S3p11, S3p21, S3p31,labels = c("A", "B","C"), ncol=3, common.legend = TRUE, font.label = list(size=18),  legend="right")
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_Sim1_qvalue_0.01vs50_initial_letters.pdf",width = 21,height = 7,pointsize = 12)
ggarrange(S3p11, S3p21, S3p31,labels = c("A", "B","C"), ncol=3, common.legend = TRUE, font.label = list(size=18),  legend="right")
dev.off()

tiff(filename = "PCA_RDA_DL_MLP_Scan_Sim1_14_7_qvalue_0.01vs50_initial_letters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
ggarrange(S3p11, S3p21, S3p31,labels = c("A", "B","C"), ncol=3, common.legend = TRUE,font.label = list(size=18),   legend="right")
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_Sim1_14_7_qvalue_0.01vs50_initial_letters.pdf",width = 14,height = 7,pointsize = 12)
ggarrange(S3p11, S3p21, S3p31,labels = c("A", "B","C"), ncol=3, common.legend = TRUE,font.label = list(size=18),   legend="right")
dev.off()

save.image(file = "SimData_27_10_2020.RData")





#######################################################################################################
#### Comparison DL_mlp PCADAPT & RDA over the 100 simulations ####
#### DL_mlp
### split the data into 100 simulations
DLdatasim=lapply(samplesim,function(x) {DLdata[,x]})
q.values_DLmlp_simus <- lapply(DLdatasim, function(x) {DLqvalues(x, 10)[,2]})
adj.pvalues_DLmlp_simus <- lapply(DLdatasim, function(x) {DLqvalues(x, 10)[,3]})
### adj.pvalue is another form of q.value

#RDA
fich_simu<-list.files(recursive = T, path="./", pattern="geno.pcadapt", full.names = TRUE)
genos<-lapply(fich_simu, function(x) t(read.table(x)))
env<-read.csv("./sim1/sim1.csv")#[,1:10]

RDA_res <- lapply(genos, function(x) {rda(x ~ env$envir1 + env$envir2 + env$envir3 + env$envir4 + env$envir5 + env$envir6 + env$envir7 + env$envir8 + env$envir9 + env$envir10,  env)})
q.values_rdadapt_simus <- lapply(RDA_res, function(x) {rdadapt(x, 10)[,2]})### the original article k =5

# PCADAPT
genos_pcadapt<-lapply(fich_simu, function(x) read.pcadapt(x, type = "pcadapt"))
pcadapt <- lapply(genos_pcadapt, function(x) pcadapt(x, K=10, method="mahalanobis",  min.maf = 0.05, ploidy = 2))
q.values_pcadapt_simus <- lapply(pcadapt, function(x) qvalue(x$pvalues)$qvalues)

save.image(file = "DL_Scan_PCA_RDA_DL_mlp_100sim.RData")
## Power = number of true positive found by the analysis for each QTL
Loci<-rep("Neutral", 1000)
Loci[c(1,11,21,31,41,51,61,71,81,91)]<-"QTL1"
Loci[c(101,111,121,131,141,151,161,171,181,191)]<-"QTL2"
Loci[c(201,211,221,231,241,251,261,271,281,291)]<-"QTL3"
### first set the cutoff as 0.01, number of loci detected by a method per sim
Tpos_pcadapt_QTL1<-unlist(lapply(q.values_pcadapt_simus, function(x) length(which(x < 0.01 & Loci=="QTL1"))))
Tpos_pcadapt_QTL2<-unlist(lapply(q.values_pcadapt_simus, function(x) length(which(x < 0.01 & Loci=="QTL2"))))
Tpos_pcadapt_QTL3<-unlist(lapply(q.values_pcadapt_simus, function(x) length(which(x < 0.01 & Loci=="QTL3"))))

Tpos_rdadapt_QTL1<-unlist(lapply(q.values_rdadapt_simus, function(x) length(which(x < 0.01 & Loci=="QTL1"))))
Tpos_rdadapt_QTL2<-unlist(lapply(q.values_rdadapt_simus, function(x) length(which(x < 0.01 & Loci=="QTL2"))))
Tpos_rdadapt_QTL3<-unlist(lapply(q.values_rdadapt_simus, function(x) length(which(x < 0.01 & Loci=="QTL3"))))

Tpos_DLmlp_QTL1<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 0.01  & Loci=="QTL1"))))
Tpos_DLmlp_QTL2<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 0.01 & Loci=="QTL2"))))
Tpos_DLmlp_QTL3<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 0.01  & Loci=="QTL3"))))

## adju pvalue, do these later to compare the q.values and then choose one form
Tpos_DLmlp_QTL1<-unlist(lapply(adj.pvalues_DLmlp_simus, function(x) length(which(x < 0.01  & Loci=="QTL1"))))
Tpos_DLmlp_QTL2<-unlist(lapply(adj.pvalues_DLmlp_simus, function(x) length(which(x < 0.01 & Loci=="QTL2"))))
Tpos_DLmlp_QTL3<-unlist(lapply(adj.pvalues_DLmlp_simus, function(x) length(which(x < 0.01  & Loci=="QTL3"))))

## False discovery rate 
Fpos_pcadapt<-unlist(lapply(q.values_pcadapt_simus, function(x) length(which(x < 0.01 & Loci=="Neutral"))))
Fpos_rdadapt<-unlist(lapply(q.values_rdadapt_simus, function(x) length(which(x < 0.01 & Loci=="Neutral"))))

Fpos_DLmlp<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x <0.01 & Loci=="Neutral"))))

Fpos_DLmlp<-unlist(lapply(adj.pvalues_DLmlp_simus, function(x) length(which(x <0.01 & Loci=="Neutral"))))

### t test to see the difference
t.test(Tpos_DLmlp_QTL1/10,Tpos_rdadapt_QTL1/10)
t.test(Tpos_DLmlp_QTL1/10,Tpos_pcadapt_QTL1/10)
t.test(Tpos_rdadapt_QTL1/10,Tpos_pcadapt_QTL1/10)


t.test(Tpos_DLmlp_QTL2/10,Tpos_rdadapt_QTL2/10)
t.test(Tpos_DLmlp_QTL2/10,Tpos_pcadapt_QTL2/10)
t.test(Tpos_rdadapt_QTL2/10,Tpos_pcadapt_QTL2/10)

t.test(Tpos_DLmlp_QTL3/10,Tpos_rdadapt_QTL3/10)
t.test(Tpos_DLmlp_QTL3/10,Tpos_pcadapt_QTL3/10)
t.test(Tpos_rdadapt_QTL3/10,Tpos_pcadapt_QTL3/10)


Allpositivespcadapt0.01=unlist(lapply(q.values_pcadapt_simus, function(x) length(which(x < 0.01))))
Allpositivesrda0.01=unlist(lapply(q.values_rdadapt_simus, function(x) length(which(x < 0.01))))

AllpositivesDLmlp0.01=unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 0.01))))

AllpositivesDLmlp0.01=unlist(lapply(adj.pvalues_DLmlp_simus, function(x) length(which(x < 0.01))))

### True Discovery rate
Tpos<-data.frame(Tpos=c(Tpos_pcadapt_QTL1/10, Tpos_pcadapt_QTL2/10, Tpos_pcadapt_QTL3/10, Tpos_rdadapt_QTL1/10, Tpos_rdadapt_QTL2/10, Tpos_rdadapt_QTL3/10,Tpos_DLmlp_QTL1/10, Tpos_DLmlp_QTL2/10, Tpos_DLmlp_QTL3/10), Method=c(rep("PCA", 300), rep("RDA", 300),rep("DL_mlp", 300)), Loci=c(rep("QTL1", 100), rep("QTL2", 100), rep("QTL3", 100), rep("QTL1", 100), rep("QTL2", 100), rep("QTL3", 100),rep("QTL1", 100), rep("QTL2", 100), rep("QTL3", 100)))
Tpos_meansd<-cbind(aggregate(Tpos$Tpos, by=list(Tpos$Method,Tpos$Loci), FUN = "mean"), sd=aggregate(Tpos$Tpos, by=list(Tpos$Method,Tpos$Loci), FUN = "sd")$x)
colnames(Tpos_meansd)[1:3]<-c("Method", "Loci", "Rate")




## Visualisation
Tposplot <- ggplot(Tpos_meansd, aes(y=Rate, x=Loci, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'lightgray')) +
  ylab("Percentage of True Positive") +
  guides(fill=FALSE) +
  xlab("") +
  theme_classic()
### with letters

Tposplot <- ggplot(Tpos_meansd, aes(y=Rate, x=Loci, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("True Discovery Rate") +
  guides(fill=FALSE) +
  xlab("") +
  theme(axis.text=element_text(size=12)) +geom_text(
    aes(label = c("a","c","b","a","c","b","a","c","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)



#### update the FDR :Number of false positives/Total number of all positives
Fpos<-data.frame(Fpos=c(Fpos_pcadapt/Allpositivespcadapt0.01, Fpos_rdadapt/Allpositivesrda0.01,Fpos_DLmlp/AllpositivesDLmlp0.01), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DL_MLP", 100)))
Fpos<- cbind(aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "mean"), sd=aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "sd")$x)
colnames(Fpos)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp/AllpositivesDLmlp0.01,Fpos_rdadapt/Allpositivesrda0.01)
t.test(Fpos_DLmlp/AllpositivesDLmlp0.01,Fpos_pcadapt/Allpositivespcadapt0.01)
t.test(Fpos_rdadapt/Allpositivesrda0.01,Fpos_pcadapt/Allpositivespcadapt0.01)


Fposplot <- ggplot(Fpos, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'lightgray')) +
  ylab("False Discovery Rate") +
  ylim(0,2) +  
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
### with letters
Fposplot <- ggplot(Fpos, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,1) +  
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("a","b","c"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5,size=6)

#### all 0.01
tiff(filename = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvalue_MLPadjp-value_all_0.01_withletters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot,Fposplot, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_aPCA_RDA_qvalue_MLPdjp-value_all_0.01_withletters.pdf",width = 14,height = 7,pointsize = 12)
plot_grid(Tposplot,Fposplot, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()

save(Fpos,Tpos, Tpos_meansd,Fpos_pcadapt,Allpositivespcadapt0.01, Fpos_rdadapt,Allpositivesrda0.01,Fpos_DLmlp,AllpositivesDLmlp0.01,Tposplot,Fposplot,file = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjp0.01_withletters.RData")



save.image(file = "Simulations_PCA_RDA_DL_mlpScan_25_june_2020.RData")
###############################################################################################
### cutoff of 1e-25/ 1e-50 1e-130, 100 no difference between fasle discovery rate, 125 MLP is significantly lower than the other two 

Tpos_pcadapt_QTL1<-unlist(lapply(q.values_pcadapt_simus, function(x) length(which(x < 1e-100 & Loci=="QTL1"))))
Tpos_pcadapt_QTL2<-unlist(lapply(q.values_pcadapt_simus, function(x) length(which(x < 1e-100 & Loci=="QTL2"))))
Tpos_pcadapt_QTL3<-unlist(lapply(q.values_pcadapt_simus, function(x) length(which(x < 1e-100 & Loci=="QTL3"))))

Tpos_rdadapt_QTL1<-unlist(lapply(q.values_rdadapt_simus, function(x) length(which(x < 1e-100 & Loci=="QTL1"))))
Tpos_rdadapt_QTL2<-unlist(lapply(q.values_rdadapt_simus, function(x) length(which(x < 1e-100& Loci=="QTL2"))))
Tpos_rdadapt_QTL3<-unlist(lapply(q.values_rdadapt_simus, function(x) length(which(x <1e-100& Loci=="QTL3"))))

Tpos_DLmlp_QTL1<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-100  & Loci=="QTL1"))))
Tpos_DLmlp_QTL2<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-100 & Loci=="QTL2"))))
Tpos_DLmlp_QTL3<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-100  & Loci=="QTL3"))))

#### skip and do this later again

Tpos_DLmlp_QTL1<-unlist(lapply(adj.pvalues_DLmlp_simus, function(x) length(which(x < 1e-100 & Loci=="QTL1"))))
Tpos_DLmlp_QTL2<-unlist(lapply(adj.pvalues_DLmlp_simus, function(x) length(which(x < 1e-100 & Loci=="QTL2"))))
Tpos_DLmlp_QTL3<-unlist(lapply(adj.pvalues_DLmlp_simus, function(x) length(which(x < 1e-100  & Loci=="QTL3"))))


##### skip and do this later again
Tpos_DLmlp_QTL1<-unlist(lapply(adj.pvalues_DLmlp_simus, function(x) length(which(x < 1e-115  & Loci=="QTL1"))))
Tpos_DLmlp_QTL2<-unlist(lapply(adj.pvalues_DLmlp_simus, function(x) length(which(x < 1e-115 & Loci=="QTL2"))))
Tpos_DLmlp_QTL3<-unlist(lapply(adj.pvalues_DLmlp_simus, function(x) length(which(x < 1e-115  & Loci=="QTL3"))))


####################26-10-2020, seting a threshold of 25, 50, which are between 10-100 to see the false discovery rate
### in oder to make use of the existing script, I complete each process for each threshold, and returned do for another threshold,or else you will need copy and changed the names and repeat 
### Below scripts were done each time and returned to repeat again, there were not done sequentially , users should understand the above meaning

Tpos_DLmlp_QTL1<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-10  & Loci=="QTL1"))))
Tpos_DLmlp_QTL2<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-10 & Loci=="QTL2"))))
Tpos_DLmlp_QTL3<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-10  & Loci=="QTL3"))))


Tpos_DLmlp_QTL1<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-20  & Loci=="QTL1"))))
Tpos_DLmlp_QTL2<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-20 & Loci=="QTL2"))))
Tpos_DLmlp_QTL3<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-20  & Loci=="QTL3"))))

Tpos_DLmlp_QTL1<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-25  & Loci=="QTL1"))))
Tpos_DLmlp_QTL2<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-25 & Loci=="QTL2"))))
Tpos_DLmlp_QTL3<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-25  & Loci=="QTL3"))))


Tpos_DLmlp_QTL1<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-30  & Loci=="QTL1"))))
Tpos_DLmlp_QTL2<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-30 & Loci=="QTL2"))))
Tpos_DLmlp_QTL3<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-30  & Loci=="QTL3"))))


Tpos_DLmlp_QTL1<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-40  & Loci=="QTL1"))))
Tpos_DLmlp_QTL2<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-40 & Loci=="QTL2"))))
Tpos_DLmlp_QTL3<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-40  & Loci=="QTL3"))))

Tpos_DLmlp_QTL1<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-50  & Loci=="QTL1"))))
Tpos_DLmlp_QTL2<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-50 & Loci=="QTL2"))))
Tpos_DLmlp_QTL3<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-50  & Loci=="QTL3"))))


## False discovery rate 
Fpos_pcadapt<-unlist(lapply(q.values_pcadapt_simus, function(x) length(which(x < 1e-100 & Loci=="Neutral"))))
Fpos_rdadapt<-unlist(lapply(q.values_rdadapt_simus, function(x) length(which(x < 1e-100 & Loci=="Neutral"))))

Fpos_DLmlp<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x <1e-130 & Loci=="Neutral"))))

Fpos_DLmlp<-unlist(lapply(adj.pvalues_DLmlp_simus, function(x) length(which(x <1e-100 & Loci=="Neutral"))))

Fpos_DLmlp<-unlist(lapply(adj.pvalues_DLmlp_simus, function(x) length(which(x <1e-115 & Loci=="Neutral"))))

###25-10-2020 remember do this each time rather than sequentially

Fpos_DLmlp<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x <1e-10 & Loci=="Neutral"))))
Fpos_DLmlp<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x <1e-20 & Loci=="Neutral"))))

Fpos_DLmlp<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x <1e-25 & Loci=="Neutral"))))

Fpos_DLmlp<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x <1e-30 & Loci=="Neutral"))))
Fpos_DLmlp<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x <1e-40 & Loci=="Neutral"))))

Fpos_DLmlp<-unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x <1e-50 & Loci=="Neutral"))))

### t test to see the difference
t.test(Tpos_DLmlp_QTL1/10,Tpos_rdadapt_QTL1/10)
t.test(Tpos_DLmlp_QTL1/10,Tpos_pcadapt_QTL1/10)
t.test(Tpos_rdadapt_QTL1/10,Tpos_pcadapt_QTL1/10)


t.test(Tpos_DLmlp_QTL2/10,Tpos_rdadapt_QTL2/10)
t.test(Tpos_DLmlp_QTL2/10,Tpos_pcadapt_QTL2/10)
t.test(Tpos_rdadapt_QTL2/10,Tpos_pcadapt_QTL2/10)

t.test(Tpos_DLmlp_QTL3/10,Tpos_rdadapt_QTL3/10)
t.test(Tpos_DLmlp_QTL3/10,Tpos_pcadapt_QTL3/10)
t.test(Tpos_rdadapt_QTL3/10,Tpos_pcadapt_QTL3/10)


Allpositivespcadapt100=unlist(lapply(q.values_pcadapt_simus, function(x) length(which(x < 1e-100))))
Allpositivesrda100=unlist(lapply(q.values_rdadapt_simus, function(x) length(which(x < 1e-100))))

AllpositivesDLmlp100=unlist(lapply(adj.pvalues_DLmlp_simus, function(x) length(which(x < 1e-100))))

AllpositivesDLmlp115=unlist(lapply(adj.pvalues_DLmlp_simus, function(x) length(which(x < 1e-115))))
##########26-10-2020
AllpositivesDLmlp10=unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-10))))
AllpositivesDLmlp20=unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-20))))
AllpositivesDLmlp25=unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-25))))

AllpositivesDLmlp30=unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-30))))
AllpositivesDLmlp40=unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-40))))
AllpositivesDLmlp50=unlist(lapply(q.values_DLmlp_simus, function(x) length(which(x < 1e-50))))

### True Discovery rate
Tpos<-data.frame(Tpos=c(Tpos_pcadapt_QTL1/10, Tpos_pcadapt_QTL2/10, Tpos_pcadapt_QTL3/10, Tpos_rdadapt_QTL1/10, Tpos_rdadapt_QTL2/10, Tpos_rdadapt_QTL3/10,Tpos_DLmlp_QTL1/10, Tpos_DLmlp_QTL2/10, Tpos_DLmlp_QTL3/10), Method=c(rep("PCA", 300), rep("RDA", 300),rep("DL_mlp", 300)), Loci=c(rep("QTL1", 100), rep("QTL2", 100), rep("QTL3", 100), rep("QTL1", 100), rep("QTL2", 100), rep("QTL3", 100),rep("QTL1", 100), rep("QTL2", 100), rep("QTL3", 100)))
Tpos_meansd<-cbind(aggregate(Tpos$Tpos, by=list(Tpos$Method,Tpos$Loci), FUN = "mean"), sd=aggregate(Tpos$Tpos, by=list(Tpos$Method,Tpos$Loci), FUN = "sd")$x)
colnames(Tpos_meansd)[1:3]<-c("Method", "Loci", "Rate")




## Visualisation
Tposplot <- ggplot(Tpos_meansd, aes(y=Rate, x=Loci, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'lightgray')) +
  ylab("Percentage of True Positive") +
  guides(fill=FALSE) +
  xlab("") +
  theme_classic()


### with letters

Tposplot <- ggplot(Tpos_meansd, aes(y=Rate, x=Loci, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("Percentage of True Positive") +
  guides(fill=FALSE) +
  xlab("") +
  theme_classic() +geom_text(
    aes(label = c("a","c","b","a","c","b","a","c","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5)

###########25-20-2020

#### pca-rda 0.01, DL_MLP 1e-10

Fpos<-data.frame(Fpos=c(Fpos_pcadapt/Allpositivespcadapt0.01, Fpos_rdadapt/Allpositivesrda0.01,Fpos_DLmlp/AllpositivesDLmlp10), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DL_MLP", 100)))
Fpos<- cbind(aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "mean"), sd=aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "sd")$x)
colnames(Fpos)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp/AllpositivesDLmlp10,Fpos_rdadapt/Allpositivesrda0.01)
t.test(Fpos_DLmlp/AllpositivesDLmlp10,Fpos_pcadapt/Allpositivespcadapt0.01)
t.test(Fpos_rdadapt/Allpositivesrda0.01,Fpos_pcadapt/Allpositivespcadapt0.01)

#### pca-rda 0.01, DL_MLP 1e-20

Fpos<-data.frame(Fpos=c(Fpos_pcadapt/Allpositivespcadapt0.01, Fpos_rdadapt/Allpositivesrda0.01,Fpos_DLmlp/AllpositivesDLmlp20), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DL_MLP", 100)))
Fpos<- cbind(aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "mean"), sd=aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "sd")$x)
colnames(Fpos)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp/AllpositivesDLmlp20,Fpos_rdadapt/Allpositivesrda0.01)
t.test(Fpos_DLmlp/AllpositivesDLmlp20,Fpos_pcadapt/Allpositivespcadapt0.01)
t.test(Fpos_rdadapt/Allpositivesrda0.01,Fpos_pcadapt/Allpositivespcadapt0.01)


#### pca-rda 0.01, DL_MLP 1e-25

Fpos<-data.frame(Fpos=c(Fpos_pcadapt/Allpositivespcadapt0.01, Fpos_rdadapt/Allpositivesrda0.01,Fpos_DLmlp/AllpositivesDLmlp25), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DL_MLP", 100)))
Fpos<- cbind(aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "mean"), sd=aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "sd")$x)
colnames(Fpos)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp/AllpositivesDLmlp25,Fpos_rdadapt/Allpositivesrda0.01)
t.test(Fpos_DLmlp/AllpositivesDLmlp25,Fpos_pcadapt/Allpositivespcadapt0.01)
t.test(Fpos_rdadapt/Allpositivesrda0.01,Fpos_pcadapt/Allpositivespcadapt0.01)
#### pca-rda 0.01, DL_MLP 1e-30

Fpos<-data.frame(Fpos=c(Fpos_pcadapt/Allpositivespcadapt0.01, Fpos_rdadapt/Allpositivesrda0.01,Fpos_DLmlp/AllpositivesDLmlp30), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DL_MLP", 100)))
Fpos<- cbind(aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "mean"), sd=aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "sd")$x)
colnames(Fpos)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp/AllpositivesDLmlp30,Fpos_rdadapt/Allpositivesrda0.01)
t.test(Fpos_DLmlp/AllpositivesDLmlp30,Fpos_pcadapt/Allpositivespcadapt0.01)
t.test(Fpos_rdadapt/Allpositivesrda0.01,Fpos_pcadapt/Allpositivespcadapt0.01)




#### pca-rda 0.01, DL_MLP 1e-40

Fpos<-data.frame(Fpos=c(Fpos_pcadapt/Allpositivespcadapt0.01, Fpos_rdadapt/Allpositivesrda0.01,Fpos_DLmlp/AllpositivesDLmlp40), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DL_MLP", 100)))
Fpos<- cbind(aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "mean"), sd=aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "sd")$x)
colnames(Fpos)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp/AllpositivesDLmlp40,Fpos_rdadapt/Allpositivesrda0.01)
t.test(Fpos_DLmlp/AllpositivesDLmlp40,Fpos_pcadapt/Allpositivespcadapt0.01)
t.test(Fpos_rdadapt/Allpositivesrda0.01,Fpos_pcadapt/Allpositivespcadapt0.01)

#### pca-rda 0.01, DL_MLP 1e-50

Fpos<-data.frame(Fpos=c(Fpos_pcadapt/Allpositivespcadapt0.01, Fpos_rdadapt/Allpositivesrda0.01,Fpos_DLmlp/AllpositivesDLmlp50), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DL_MLP", 100)))
Fpos<- cbind(aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "mean"), sd=aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "sd")$x)
colnames(Fpos)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp/AllpositivesDLmlp50,Fpos_rdadapt/Allpositivesrda0.01)
t.test(Fpos_DLmlp/AllpositivesDLmlp50,Fpos_pcadapt/Allpositivespcadapt0.01)
t.test(Fpos_rdadapt/Allpositivesrda0.01,Fpos_pcadapt/Allpositivespcadapt0.01)

#### pca-rda 0.01, DL_MLP 1e-100

Fpos<-data.frame(Fpos=c(Fpos_pcadapt/Allpositivespcadapt0.01, Fpos_rdadapt/Allpositivesrda0.01,Fpos_DLmlp/AllpositivesDLmlp100), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DL_MLP", 100)))
Fpos<- cbind(aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "mean"), sd=aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "sd")$x)
colnames(Fpos)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp/AllpositivesDLmlp100,Fpos_rdadapt/Allpositivesrda0.01)
t.test(Fpos_DLmlp/AllpositivesDLmlp100,Fpos_pcadapt/Allpositivespcadapt0.01)
t.test(Fpos_rdadapt/Allpositivesrda0.01,Fpos_pcadapt/Allpositivespcadapt0.01)

#### pca-rda 0.01, DL_MLP 1e-125

Fpos<-data.frame(Fpos=c(Fpos_pcadapt/Allpositivespcadapt0.01, Fpos_rdadapt/Allpositivesrda0.01,Fpos_DLmlp/AllpositivesDLmlp115), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DL_MLP", 100)))
Fpos<- cbind(aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "mean"), sd=aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "sd")$x)
colnames(Fpos)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp/AllpositivesDLmlp115,Fpos_rdadapt/Allpositivesrda0.01)
t.test(Fpos_DLmlp/AllpositivesDLmlp115,Fpos_pcadapt/Allpositivespcadapt0.01)
t.test(Fpos_rdadapt/Allpositivesrda0.01,Fpos_pcadapt/Allpositivespcadapt0.01)



###### all 1e-100
Fpos<-data.frame(Fpos=c(Fpos_pcadapt/Allpositivespcadapt100, Fpos_rdadapt/Allpositivesrda100,Fpos_DLmlp/AllpositivesDLmlp100), Method=c(rep("PCA", 100), rep("RDA", 100),rep("DL_MLP", 100)))
Fpos<- cbind(aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "mean"), sd=aggregate(Fpos$Fpos, by=list(Fpos$Method), FUN = "sd")$x)
colnames(Fpos)[1:2]<-c("Method", "Rate")

t.test(Fpos_DLmlp/AllpositivesDLmlp100,Fpos_rdadapt/Allpositivesrda100)
t.test(Fpos_DLmlp/AllpositivesDLmlp100,Fpos_pcadapt/Allpositivespcadapt100)
t.test(Fpos_rdadapt/Allpositivesrda100,Fpos_pcadapt/Allpositivespcadapt100)






Fposplot <- ggplot(Fpos, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'lightgray')) +
  ylab("False Discovery Rate") +
  ylim(0,0.1) +  
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

Fposplot <- ggplot(Fpos, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'lightgray')) +
  ylab("False Discovery Rate") +
  ylim(0,2) +  
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
### with letters 0.01, q-vlue =-10
Fposplot <- ggplot(Fpos, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,1) + 
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("a","b","c"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5)
### with letters 0.01, q-vlue =-20
Fposplot <- ggplot(Fpos, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,1) + 
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("a","b","c"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5)
### with letters 0.01, q-vlue =-30
Fposplot <- ggplot(Fpos, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,1) + 
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("a","b","c"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5)
### with letters 0.01, q-vlue =-40
Fposplot <- ggplot(Fpos, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,1) + 
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("a","b","c"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5)
### with letters 0.01, q-vlue =-50
Fposplot <- ggplot(Fpos, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,1) + 
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("a","b","c"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5)
### with letters 0.01, q-vlue =-50
Fposplot <- ggplot(Fpos, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,1) + 
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("a","a","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5)
### with letters 0.01, q-vlue =-100
Fposplot <- ggplot(Fpos, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,0.2) + 
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("b","a","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5)


### with letters q-vlue =-115
Fposplot <- ggplot(Fpos, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,0.2) +  
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +geom_text(
    aes(label = c("c","a","b"),y=Rate+sd),
    position = position_dodge(0.9),
    vjust = -0.5)


save(Fpos,Tpos, Fpos_pcadapt,Allpositivespcadapt0.01, Fpos_rdadapt,Allpositivesrda0.01,Fpos_DLmlp,AllpositivesDLmlp10,Tposplot,Fposplot,file = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe10_withletters.RData")
save(Fpos,Tpos, Fpos_pcadapt,Allpositivespcadapt0.01, Fpos_rdadapt,Allpositivesrda0.01,Fpos_DLmlp,AllpositivesDLmlp20,Tposplot,Fposplot,file = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe20_withletters.RData")
save(Fpos,Tpos, Fpos_pcadapt,Allpositivespcadapt0.01, Fpos_rdadapt,Allpositivesrda0.01,Fpos_DLmlp,AllpositivesDLmlp25,Tposplot,Fposplot,file = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe25_withletters.RData")

save(Fpos,Tpos, Fpos_pcadapt,Allpositivespcadapt0.01, Fpos_rdadapt,Allpositivesrda0.01,Fpos_DLmlp,AllpositivesDLmlp30,Tposplot,Fposplot,file = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe30_withletters.RData")
save(Fpos,Tpos, Fpos_pcadapt,Allpositivespcadapt0.01, Fpos_rdadapt,Allpositivesrda0.01,Fpos_DLmlp,AllpositivesDLmlp40,Tposplot,Fposplot,file = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe40_withletters.RData")


save(Fpos,Tpos, Fpos_pcadapt,Allpositivespcadapt0.01, Fpos_rdadapt,Allpositivesrda0.01,Fpos_DLmlp,AllpositivesDLmlp100,Tposplot,Fposplot,file = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe100_withletters.RData")

save(Fpos,Tpos, Fpos_pcadapt,Allpositivespcadapt0.01, Fpos_rdadapt,Allpositivesrda0.01,Fpos_DLmlp,AllpositivesDLmlp115,Tposplot,Fposplot,file = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe115_withletters.RData")

###### with letters all q-vlue =-100
### with letters

Tposplot <- ggplot(Tpos_meansd, aes(y=Rate, x=Loci, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("Percentage of True Positive") +
  guides(fill=FALSE) +
  xlab("") +
  theme_classic() 


Fposplot <- ggplot(Fpos, aes(y=Rate, x=1, fill=Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Rate-sd, ymax=Rate+sd), width=.1, position = position_dodge(.9)) +
  scale_fill_manual(values=c("blue",'black', 'gray')) +
  ylab("False Discovery Rate") +
  ylim(0,1) + 
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())


save(Fpos,Tpos, Fpos_pcadapt,Allpositivespcadapt100, Fpos_rdadapt,Allpositivesrda100,Fpos_DLmlp,AllpositivesDLmlp100,Tposplot,Fposplot,file = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvalee100MLP_adjpe100withletters.RData")


###  0.01. vs -10
tiff(filename = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe10_withletters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot,Fposplot, labels = c("a","b"), align = "h", rel_widths = c(1, .7))
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe10_withletters.pdf",width = 14,height = 7,pointsize = 12)
plot_grid(Tposplot,Fposplot, labels = c("a","b"), align = "h", rel_widths = c(1, .7))
dev.off()

###  0.01. vs -20
tiff(filename = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe20_withletters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot,Fposplot, labels = c("a","b"), align = "h", rel_widths = c(1, .7))
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe20_withletters.pdf",width = 14,height = 7,pointsize = 12)
plot_grid(Tposplot,Fposplot, labels = c("a","b"), align = "h", rel_widths = c(1, .7))
dev.off()
###  0.01. vs -30
tiff(filename = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe30_withletters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot,Fposplot, labels = c("a","b"), align = "h", rel_widths = c(1, .7))
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe30_withletters.pdf",width = 14,height = 7,pointsize = 12)
plot_grid(Tposplot,Fposplot, labels = c("a","b"), align = "h", rel_widths = c(1, .7))
dev.off()

###  0.01. vs -40
tiff(filename = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe40_withletters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot,Fposplot, labels = c("a","b"), align = "h", rel_widths = c(1, .7))
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe40_withletters.pdf",width = 14,height = 7,pointsize = 12)
plot_grid(Tposplot,Fposplot, labels = c("a","b"), align = "h", rel_widths = c(1, .7))
dev.off()

###  0.01. vs -50
tiff(filename = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe50_withletters_1.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot,Fposplot, labels = c("a","b"), align = "h", rel_widths = c(1, .7))
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe50_withletters_1.pdf",width = 14,height = 7,pointsize = 12)
plot_grid(Tposplot,Fposplot, labels = c("a","b"), align = "h", rel_widths = c(1, .7))
dev.off()

save(Fpos,Tpos, Fpos_pcadapt,Allpositivespcadapt0.01, Fpos_rdadapt,Allpositivesrda0.01,Fpos_DLmlp,AllpositivesDLmlp0.01,Tposplot,Fposplot,file = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvalee0.01MLP_adjpe50withletters.RData")

### pvalue all -100
tiff(filename = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale100MLP_adjpe100_withletters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot,Fposplot, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale100MLP_adjpe100_withletters.pdf",width = 14,height = 7,pointsize = 12)
plot_grid(Tposplot,Fposplot, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()

###  0.01. vs -100
tiff(filename = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe100_withletters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot,Fposplot, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe100_withletters.pdf",width = 14,height = 7,pointsize = 12)
plot_grid(Tposplot,Fposplot, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()

###  0.01. vs -115
tiff(filename = "PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe115_withletters.tif",width = 14,height = 7,units = "in",pointsize = 12,res=600,compression = "zip+p")
plot_grid(Tposplot,Fposplot, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()
pdf("PCA_RDA_DL_MLP_Scan_True_and False_discovery_rats_100simus_PCA_RDA_qvale0.01MLP_adjpe115_withletters.pdf",width = 14,height = 7,pointsize = 12)
plot_grid(Tposplot,Fposplot, labels = "AUTO", align = "h", rel_widths = c(1, .7))
dev.off()

save.image(file = "Simulations_PCA_RDA_DL_mlpScan_10_Agu_2020.RData")




