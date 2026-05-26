#cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC/fMRI_coupling/PopulationAnalysis;module load r/4.4.0;R
iii=53
library(pls)
library(data.table)
fileT=Sys.glob('G360_PC/ph12/PCs*.txt')
A=fread(fileT[iii],data.table=F)
A1=fread(gsub('/ph12/','/ph3/',fileT[iii]),data.table=F)
A=rbind(A,A1)
rm(list=c('A1'))
outfile=gsub('_100.txt','',gsub('G360_PC/ph12/PCs_','',fileT[iii]))
rownames(A)=A[,1];A=A[,-1]
pp=dim(A)[2]
for(ii in 1:dim(A)[2])
{
print(dim(A)[2]-ii)
med0=median(A[,ii],na.rm=T)
mad0=mad(A[,ii],na.rm=T)
A[which(abs(A[,ii]-med0)>5*mad0),ii]=NA
indtemp=which(is.na(A[,ii]))
if(length(indtemp)>0)A[indtemp,ii]=median(A[,ii],na.rm=T)
}
B_train=fread('../../../PCA_ELReg/FA_UKB_Phase12/Mean.txt',data.table=F)
rownames(B_train)=B_train[,1];B_train=B_train[,-1]
B_test=fread('../../../PCA_ELReg/FA_UKB_Phase3_proj/Mean.txt',data.table=F)
rownames(B_test)=B_test[,1];B_test=B_test[,-1]
UID_train=intersect(rownames(A),rownames(B_train))
UID_test=intersect(rownames(A),rownames(B_test))
A_train=A[UID_train,];B_train=B_train[UID_train,]
A_test=A[UID_test,];B_test=B_test[UID_test,]
#attr(temp1,"scaled:center");attr(temp1,"scaled:scale")
A_train <- scale(A_train, center = TRUE, scale = TRUE)
B_train <- scale(B_train, center = TRUE, scale = TRUE)
A_test  <- scale(A_test, center = attr(A_train, "scaled:center"), scale = attr(A_train, "scaled:scale"))
B_test  <- scale(B_test, center = attr(B_train, "scaled:center"), scale = attr(B_train, "scaled:scale"))
M <- rep(NA, pp);names(M)=colnames(A)
names(M) <- colnames(A)
ncomp0=4
A_pred <- A_test;A_pred[T]=NA
for (ii in 1:pp) {
  print(pp-ii)
  pls_model <- plsr(A_train[,ii] ~ B_train, ncomp = ncomp0, validation = "none")  # choose ncomp by CV if needed
  A_pred[,ii] <- predict(pls_model, newdata = B_test, ncomp = ncomp0)[,1,1]
  M[ii] <- cor(A_pred[,ii], A_test[, ii])
  print(M[ii])
}
write.table(M,file=paste0('G360_PC_out/',outfile,'_cor.txt'),quote=F)
A_pred1=t(t(A_pred)*attr(A_train,"scaled:scale")+attr(A_train,"scaled:center"));colnames(A_pred1)=colnames(A)
write.table(A_pred1,file=paste0('G360_PC_out/',outfile,'_pred.txt'),quote=F)
pcnum=read.table(paste0('G360_PC/ph12/PC_',outfile,'.txt'))[[1]]
ind00=which(pcnum>0.9);qqmax=100;if(length(ind00)>0)qqmax=ind00[1]
qq=1:qqmax
A_pred1=A_pred1[,qq]
Eigv=read.table(paste0('G360_PC/ph12/Eigv_',outfile,'_100.txt'))
sd1=read.table(paste0('G360_PC/ph12/sd_',outfile,'.txt'))[[1]]
mu1=read.table(paste0('G360_PC/ph12/mu0_',outfile,'.txt'))[[1]]
AT=fread(paste0('G360/',outfile,'.txt.gz'),data.table=F)
rownames(AT)=AT[,1];AT=AT[,-1]
pp=dim(AT)[2]
for(ii in 1:dim(AT)[2])
{
med0=median(AT[,ii],na.rm=T)
mad0=mad(AT[,ii],na.rm=T)
AT[which(abs(AT[,ii]-med0)>5*mad0),ii]=NA
indtemp=which(is.na(AT[,ii]))
if(length(indtemp)>0)AT[indtemp,ii]=median(AT[,ii],na.rm=T)
}
A_test1=AT[rownames(A_test),]
A_train1=AT[rownames(A_train),]
A_pred1=as.matrix(A_pred1)%*%t(as.matrix(Eigv)[,qq])%*%diag(sd1^2)
A_pred1=t(t(A_pred1)+mu1)
#cor.test(A_test1[,1],A_pred1[,1])
#
cb_metrics <- function(A_pred1, A_test1, A_train1 = NULL) {
  A_pred1 <- as.matrix(A_pred1) #A_pred
  A_test1 <- as.matrix(A_test1)
  n <- nrow(A_test1)
  row_corr <- function(A, B) {
    A <- scale(A, center = TRUE, scale = FALSE)
    B <- scale(B, center = TRUE, scale = FALSE)
    num <- rowSums(A * B)
    den <- sqrt(rowSums(A^2) * rowSums(B^2))
    num / den
  }
  mu_tr <- colMeans(A_train1)
  r_rows_demean <- row_corr(sweep(A_pred1, 2, mu_tr),sweep(A_test1, 2, mu_tr))
  avgcorr_demean <- mean(r_rows_demean, na.rm = TRUE)#0.157
  #3. avgrank
  Xn <- t(scale(t(A_test1)))
  Yn <- t(scale(t(A_pred1)))
  C <- Xn %*% t(Yn)/dim(Xn)[2]  # n×n correlation matrix
  mm=rep(NA,dim(Xn)[1]);for(ii in 1:dim(Xn)[1])mm[ii]= mean(C[ii,]<C[ii,ii])
  avgrank <- mean(mm) #0.654
  avgcorr <- mean(diag(C))
  data.frame(avgcorr, avgcorr_demean, avgrank)
}
cb_metrics0=cb_metrics(A_pred1, A_test1, A_train1)
write.table(cb_metrics0,file=paste0('G360_out/',outfile,'_metrics.txt'),quote=F)
write.table(A_pred1,file=paste0('G360_out/',outfile,'_pred.txt'),quote=F)