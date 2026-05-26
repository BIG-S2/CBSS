#cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC/Aging;module load r/4.4.0;R
iii=100# 1:107
T0=100
library(pls)
library(data.table)
dimP=c(10,65,21,11)
dimP1=c(0,cumsum(dimP))
jj=which(iii<=dimP1)[1]-1
ii=iii-dimP1[jj]
PCAN=3
B_traintemp=fread('../../PCA_ELReg/FA_UKB_Phase12/test_PCs_1.txt',data.table=F)
rownames(B_traintemp)=B_traintemp[,1];B_traintemp=B_traintemp[,-1]
B_train=matrix(NA,dim(B_traintemp)[1],6090*PCAN+6090);rownames(B_train)=rownames(B_traintemp);colnames(B_train)=c(rep(NA,6090*PCAN),paste0('Mean_',1:6090))
for(ii1 in 1:6090){
if(ii1%%100==0)print(6090-ii1)
B_traintemp=fread(paste0('../../PCA_ELReg/FA_UKB_Phase12/test_PCs_',ii1,'.txt'),data.table=F)
rownames(B_traintemp)=B_traintemp[,1];B_traintemp=B_traintemp[,-1]
B_train[,((ii1-1)*PCAN+1):(ii1*PCAN)]=as.matrix(B_traintemp[,1:PCAN])
colnames(B_train)[((ii1-1)*PCAN+1):(ii1*PCAN)]=paste0('F_',ii1,'_PC',1:PCAN)
}
B_traintemp=fread('../../PCA_ELReg/FA_UKB_Phase12/Mean.txt',data.table=F)
rownames(B_traintemp)=B_traintemp[,1];B_traintemp=B_traintemp[,-1]
B_train[,(6090*PCAN+1):(6090*(PCAN+1))]=as.matrix(B_traintemp)
rm(B_traintemp)

B_testtemp=fread('../../PCA_ELReg/FA_UKB_Phase3_proj/Mean.txt',data.table=F)
rownames(B_testtemp)=B_testtemp[,1];B_testtemp=B_testtemp[,-1]
B_test=matrix(NA,dim(B_testtemp)[1],6090*PCAN+6090);rownames(B_test)=rownames(B_testtemp);colnames(B_test)=c(rep(NA,6090*PCAN),paste0('Mean_',1:6090))
for(ii1 in 1:6090){
if(ii1%%100==0)print(6090-ii1)
B_testtemp=fread(paste0('../../PCA_ELReg/FA_UKB_Phase3_proj/test_PCs_',ii1,'.txt'),data.table=F)
rownames(B_testtemp)=B_testtemp[,1];B_testtemp=B_testtemp[,-1]
B_test[,((ii1-1)*PCAN+1):(ii1*PCAN)]=as.matrix(B_testtemp[,1:PCAN])
colnames(B_test)[((ii1-1)*PCAN+1):(ii1*PCAN)]=paste0('F_',ii1,'_PC',1:PCAN)
}
B_testtemp=fread('../../PCA_ELReg/FA_UKB_Phase3_proj/Mean.txt',data.table=F)
rownames(B_testtemp)=B_testtemp[,1];B_testtemp=B_testtemp[,-1]
B_test[,(6090*PCAN+1):(6090*(PCAN+1))]=as.matrix(B_testtemp)
rm(B_testtemp)


B_train <- scale(B_train, center = TRUE, scale = TRUE)
B_test  <- scale(B_test, center = attr(B_train, "scaled:center"), scale = attr(B_train, "scaled:scale"))
#
cor_out='out/cor_pred_60900_withMean/'
system(paste0('mkdir -p ',cor_out))
fileT=Sys.glob('variables/*.csv')
fileT0=gsub('variables/','',fileT)
A=fread(fileT[jj],data.table=F)
rownames(A)=A[,1];A=A[,-1]
pp=dim(A)[2]
ncomp0=ifelse(grepl("Age|Sex|age|sex|Birth|birth", colnames(A)[ii]),15,4)
cor0=rep(NA,T0+1);
print(c(length(fileT)-jj,pp-ii))
B_train1=B_train;
B_test1=B_test
UID_train1=intersect(rownames(A)[which(!is.na(A[,ii]))],rownames(B_train1))
UID_test1=intersect(rownames(A)[which(!is.na(A[,ii]))],rownames(B_test1))
A_train1=A[UID_train1,]
A_test1=A[UID_test1,]
mu0=mean(A_train1[,ii]);sd0=sd(A_train1[,ii])
if(!grepl("Sex|sex", colnames(A)[ii])){
A_train1[,ii]=(A_train1[,ii]-mu0)/sd0
A_test1[,ii]=(A_test1[,ii]-mu0)/sd0
}
B_train1=B_train1[UID_train1,]
B_test1=B_test1[UID_test1,]
pls_temp <- plsr(A_train1[,ii] ~ B_train1, ncomp = ncomp0, validation = "none")
pred_temp <- predict(pls_temp, newdata = B_test1, ncomp = ncomp0)[,1,1]
write.table(pred_temp,file=paste0(cor_out,'/Pred_',gsub('.csv','',fileT0[jj]),'_',colnames(A)[ii],'.csv'),sep=',',quote=F,row.names=F)
corT=rep(0,dim(B_test1)[2]);for(ii1 in 1:dim(B_test1)[2])corT[ii1]=cor(pred_temp,B_test1[,ii1])
cor0[1]=cor(A_test1[,ii],pred_temp)
if(grepl("Sex|sex", colnames(A)[ii]))
{
datasex=data.frame(list(pred=as.numeric(pred_temp>0.5),true=A_test1[,ii]))
tab0=table(datasex)
cor0[1]=mean(diag(tab0)/apply(tab0,2,sum))
}
print(cor0[1])#(3PC,15):0.948;(3PC,20):0.942;(3PC,10):0.946;(5PC,15):0.951;(4PC,10):0.949; (4PC,15):0.950 (1PC,15):0.924 (2PC,15):0.940 #age 0.782
set.seed(2025)
print('Boots:')
for(ii1 in 1:T0)
{
	ind12=sample(1:dim(B_test1)[1],dim(B_test1)[1],replace=T)
	cor0[ii1+1]=cor(A_test1[ind12,ii],pred_temp[ind12])
	if(grepl("Sex|sex", colnames(A)[ii]))
	{
	datasex=data.frame(list(pred=as.numeric(pred_temp[ind12]>0.5),true=A_test1[ind12,ii]))
	tab0=table(datasex)
	cor0[ii1+1]=mean(diag(tab0)/apply(tab0,2,sum))
	}
}
write.table(corT,file=paste0(cor_out,'/Weight_',gsub('.csv','',fileT0[jj]),'_',colnames(A)[ii],'.csv'),sep=',',quote=F,row.names=F)
write.table(cor0,file=paste0(cor_out,'/Boots_',gsub('.csv','',fileT0[jj]),'_',colnames(A)[ii],'.csv'),sep=',',quote=F,row.names=F)
