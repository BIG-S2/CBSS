###################################################################################################
#Proj
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality
module load r/4.4.0;R
Template='PCA_ELReg/FA_UKB_Phase12/'
library(R.matlab)
#for(Path in c('step2_FA_UKB_Phase3','step2_FA_UKB_Retest','step2_FA_UKB_Phase4'))
#for(Path in c('step2_FA_UKB_Phase3'))
#for(Path in c('step2_FA_ABCD','step2_FA_HCP'))
#for(Path in c('step2_FA_ADNIGO2','step2_FA_HCPA','step2_FA_PING','step2_FA_PNC','step2_FA_HCPD'))
for(Path in c('step2_FA_PING'))
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
		out=paste0('PCA_ELReg/',gsub('step2_','',Path),'_proj')
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
#Mean&Max
###################################################################################################
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/
module load r/4.4.0;R
library(R.matlab)
MeanFun<-function(Path0,fun0){
    Path=paste0('step2_FA_',Path0,'/')
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
	for(ii in 2:609)
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
	out=paste0('PCA_ELReg/',gsub('\\/','',gsub('step2_','',Path)),'_proj')
	write.table(round(M,3),file=paste0(out,'/',fun0,'.txt'),quote=F,col.names=F)
}
MeanFun('ADNIGO2','mean')
MeanFun('ADNIGO2','max')
MeanFun('HCP','mean')
MeanFun('HCP','max')
MeanFun('HCPA','mean')
MeanFun('HCPA','max')
MeanFun('HCPD','mean')
MeanFun('HCPD','max')
MeanFun('PING','mean')
MeanFun('PING','max')
MeanFun('PNC','mean')
MeanFun('PNC','max')




MeanFun('ABCD','mean')
MeanFun('ABCD','max')
############################
