#cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality;module load r/4.4.0;R
library(R.matlab)
for(Path in c('step2_final/step2_FA_UKB_Phase12/','step2_final/step2_FA_UKB_Phase3/'))
{
	for(ii in 1:609)
	{
		print(609-ii)
		file=paste0(Path,ii,'_Aligned.mat')
		A=readMat(file)[[1]]#200,19204,10
		A=aperm(A,c(2,3,1))
		p=6090
		NN=10
		Vec=((ii-1)*NN+1):min(ii*NN,p)
		n=dim(A)[1]
		out=paste0('PCA_ELReg/',gsub('step2_','',Path))
		system(paste0('mkdir -p ',out))
		sign0<-function(vec){
		thr=max(abs(vec))
		ss=sign(vec)[abs(vec)>=thr]
		return(sum(ss==1)/length(ss))
		}
		p=length(Vec)
		S=10
		eigvT=array(NA,c(p,100,S))
		signT=matrix(NA,p,S)
		SVDT=MuT=matrix(NA,p,100)
		n=dim(A)[1]
		for(j in 1:p)
		{
			M=A[,j,]
			med=apply(M,2,median)
			mad=apply(M,2,mad)
			temp1=abs(M-rep(1,n)%*%t(med))/(rep(1,n)%*%t(mad))
			M[temp1>5]=NA
			mu0=apply(M,2,mean,na.rm=T)
			M[is.na(M)]=(rep(1,n)%*%t(mu0))[is.na(M)]
			L=dim(M)
			M=M-rep(1,L[1])%*%t(mu0)
			svd0=svd(M)
			eig0=svd0$d^2
			eigv=svd0$v
			SVDT[j,]=eig0
			MuT[j,]=mu0
			sign1=apply(eigv[,1:S],2,sign0)
			signT[j,]=sign1
			sign1=sign1+1e-6
			eigv[,1:S]=eigv[,1:S]*(rep(1,L[2])%*%t(sign(sign1-0.5)))
			eigvT[j,,]=eigv[,1:S]
		}
		write.table(SVDT,file=paste0(out,'/',ii,'_Eig.txt'),row.names=F,col.names=F,quote=F)
		write.table(MuT,file=paste0(out,'/',ii,'_Mu0.txt'),row.names=F,col.names=F,quote=F)
		saveRDS(eigvT,paste0(out,'/Eigv_',ii,'.rds'))
	}
}
###################################################################################################
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality
module load r/4.4.0;R
library(R.matlab)
for(Path in c('step2_final/step2_FA_UKB_Phase12/','step2_final/step2_FA_UKB_Phase3/')){
	for(ii in 1:609)
	{
		print(609-ii)
		file=paste0(Path,ii,'_Aligned.mat')
		A=readMat(file)[[1]]#200,19204,10
		A=aperm(A,c(2,3,1))
		p=6090
		NN=10
		Vec=((ii-1)*NN+1):min(ii*NN,p)
		n=dim(A)[1]
		out=paste0('PCA_ELReg/',gsub('step2_','',Path))
		muT=read.table(paste0(out,'/',ii,'_Mu0.txt'))
		system(paste0('mkdir -p ',out))
		S=10
		eigvT=readRDS(paste0(out,'/Eigv_',ii,'.rds'))
		ID=read.table(Sys.glob(paste0(Path,'*.txt')))[[1]]
		for(j in 1:length(Vec))
		{
			M=A[,j,]
			mu0=as.numeric(muT[j,])
			L=dim(M)
			M=M-rep(1,L[1])%*%t(mu0)
			M1=M%*%as.matrix(eigvT[j,,])
			rownames(M1)=ID
			write.table(M1,paste0(out,'/test_PCs_',Vec[j],'.txt'),quote=F,col.names=F)
		}
	}
}
###################################################################################################
#Proj
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality
module load r/4.4.0;R
Template='PCA_ELReg/FA_UKB_Phase12/'
library(R.matlab)
#for(Path in c('step2_final/step2_FA_UKB_Phase3','step2_final/step2_FA_UKB_Retest','step2_final/step2_FA_UKB_Phase4'))
for(Path in c('step2_final/step2_FA_UKB_Phase4'))
#for(Path in c('step2_FA_UKB_Phase3'))
#for(Path in c('step2_FA_UKB_Retest'))
{
	for(ii in 1:609)
	{
		if(ii%%10==0)print(c(609-ii,Path))
		file=paste0(Path,'/',ii,'_Aligned.mat')
		A=readMat(file)[[1]]#200,19204,10
		A=aperm(A,c(2,3,1))
		p=6090
		NN=10
		Vec=((ii-1)*NN+1):min(ii*NN,p)
		n=dim(A)[1]
		out=paste0('PCA_ELReg/',gsub('step2_final/step2_','',Path),'_proj')
		muT=read.table(paste0(Template,'/',ii,'_Mu0.txt'))
		system(paste0('mkdir -p ',out))
		S=10
		eigvT=readRDS(paste0(Template,'/Eigv_',ii,'.rds'))
		ID=read.table(Sys.glob(paste0(Path,'/*.txt')))[[1]]
		for(j in 1:length(Vec))
		{
			M=A[,j,]
			mu0=as.numeric(muT[j,])
			L=dim(M)
			M=M-rep(1,L[1])%*%t(mu0)
			M1=M%*%as.matrix(eigvT[j,,])
			rownames(M1)=ID
			write.table(M1,paste0(out,'/test_PCs_',Vec[j],'.txt'),quote=F,col.names=F)
		}
	}
}
###################################################################################################
#Reproducibility
###################################################################################################
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality
module load r/4.1.0;R
Path1='PCA_ELReg/final/FA_UKB_Phase3_proj/test_PCs_'
Path0='PCA_ELReg/final/FA_UKB_Phase3/test_PCs_'
Reprod=matrix(NA,6090,10)
for(ii in 1:6090)
{
if(ii%%100==0)print(6090-ii)
temp0=read.table(paste0(Path0,ii,'.txt'))
temp1=read.table(paste0(Path1,ii,'.txt'))
tempp=abs(cor(cbind(temp0[,-1],temp1[,-1])))
diag(tempp)=0
Reprod[ii,]=round(apply(tempp,1,max)[1:10],3)
}
write.table(Reprod,'PCA_ELReg/final/FA_UKB_Phase3_proj/Reprod_score.txt',quote=F,row.names=F,col.names=F)
######################################################################################################################################################################################################
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality
module load r/4.1.0;R
Path1='PCA_ELReg/FA_UKB_Phase12/test_PCs_'
Path0='PCA_ELReg/FA_UKB_Retest_proj/test_PCs_'
Reprod=matrix(NA,6090,10)
for(ii in 1:6090)
{
if(ii%%100==0)print(6090-ii)
temp0=read.table(paste0(Path0,ii,'.txt'))
temp1=read.table(paste0(Path1,ii,'.txt'))
rownames(temp0)=temp0[,1];rownames(temp1)=temp1[,1]
temp0=temp0[,-1];temp1=temp1[,-1]
med=apply(temp0,2,median)
mad=apply(temp0,2,mad)
ntemp=dim(temp0)[1]
temptemp1=abs(temp0-rep(1,ntemp)%*%t(med))/(rep(1,ntemp)%*%t(mad))
temp0[temptemp1>5]=NA
mu0=apply(temp0,2,mean,na.rm=T)
temp0[is.na(temp0)]=(rep(1,ntemp)%*%t(mu0))[is.na(temp0)]
med=apply(temp1,2,median)
mad=apply(temp1,2,mad)
ntemp=dim(temp1)[1]
temptemp1=abs(temp1-rep(1,ntemp)%*%t(med))/(rep(1,ntemp)%*%t(mad))
temp1[temptemp1>5]=NA
mu0=apply(temp1,2,mean,na.rm=T)
temp1[is.na(temp1)]=(rep(1,ntemp)%*%t(mu0))[is.na(temp1)]
UID=intersect(rownames(temp0),rownames(temp1))
temp0=temp0[UID,];temp1=temp1[UID,]
Reprod[ii,]=round(abs(diag(cor(cbind(temp0,temp1))[1:10,11:20])),3)
}
write.table(Reprod,'PCA_ELReg/FA_UKB_Phase3_proj/Reprod_score_test_retest.txt',quote=F,row.names=F,col.names=F)
A=read.table('PCA_ELReg/FA_UKB_Phase3_proj/Reprod_score_test_retest.txt')
indd=sort(which((A[,1]<0.6)|(A[,2]<0.6)|(A[,3]<0.6)))
tempp=A[indd,1:3]
rownames(tempp)=indd
###################################################################################################
#Contribution
###################################################################################################
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/PCA_ELReg/final/FA_UKB_Phase12 
module load r/4.4.0;R
Contri=matrix(NA,6090,10)
p=6090
NN=10
for(ii in 1:609)
{
if(ii%%100==0)print(6090-ii)
Vec=((ii-1)*NN+1):min(ii*NN,p)
temp0=read.table(paste0(ii,'_Eig.txt'))
temp00=as.matrix((t(apply(temp0,1,cumsum))/apply(temp0,1,sum))[,1:10])
Contri[Vec,]=temp00
}
Contri=round(Contri,3)
write.table(Contri,'../FA_UKB_Phase3_proj/Contri.txt',quote=F,row.names=F,col.names=F)


###################################################################################################
#Mean&Max
###################################################################################################
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/
module load r/4.4.0;R
library(R.matlab)
MeanFun<-function(Path0,fun0){
    Path=paste0('step2_final/step2_FA_',Path0,'/')
	fun1=get(fun0)
	fun0=paste0(toupper(substring(fun0, 1, 1)), substring(fun0, 2))
	p=6090
	NN=10
	ii=1
	file=paste0(Path,ii,'_Aligned.mat')
	A=readMat(file)[[1]]#200,19204,10
	A=aperm(A,c(2,3,1))
	Vec=((ii-1)*NN+1):min(ii*NN,p)
	n=dim(A)[1]
	M=matrix(NA,n,p)
	Vec=((ii-1)*NN+1):min(ii*NN,p)
	temp=apply(A,c(1,2),fun1)
	M[,Vec]=temp
	for(ii in 2:609)#609
	{
		print(609-ii)
		file=paste0(Path,ii,'_Aligned.mat')
		A=readMat(file)[[1]]#200,19204,10
		A=aperm(A,c(2,3,1))
		Vec=((ii-1)*NN+1):min(ii*NN,p)
		temp=apply(A,c(1,2),fun1)
		M[,Vec]=temp[,1:length(Vec)]
	}
	ID=read.table(Sys.glob(paste0(Path,'*.txt')))[[1]]
	rownames(M)=ID
	if(Path0=='UKB_Phase12')out=paste0('PCA_ELReg/',gsub('step2_','',Path))
	if(Path0!='UKB_Phase12')out=gsub('/_proj','_proj',paste0('PCA_ELReg/',gsub('step2_final/step2_','',Path),'_proj'))
	print(paste0(out,'/',fun0,'.txt'))
	write.table(round(M,3),file=paste0(out,'/',fun0,'.txt'),quote=F,col.names=F)
}
MeanFun('UKB_Phase12','mean')
MeanFun('UKB_Phase12','max')
MeanFun('UKB_Phase3','mean')
MeanFun('UKB_Phase3','max')
MeanFun('UKB_Phase4','mean')
MeanFun('UKB_Phase4','max')
MeanFun('UKB_Retest','mean')
MeanFun('UKB_Retest','max')
############################
