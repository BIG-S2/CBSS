#cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC/fMRI_coupling/PopulationAnalysis/;module load r/4.4.0;R
iii=2 #1:61
library(data.table)
subT=Sys.glob(paste0('G360_out/*_pred.txt'))
patT=gsub('G360_out/','',gsub('_pred.txt','',subT))
L=length(patT)
temp1=fread(subT[1],data.table=F)
rownames(temp1)=temp1[,1];temp1=temp1[,-1]
##read in the true data
INDD=((iii-1)*300+1):min(iii*300,dim(temp1)[1])
L0=length(INDD)
M1=matrix(NA,L0,64620)
#Mu=rep(NA,64620) #####
rownames(M1)=rownames(temp1)[INDD];colnames(M1)=1:64620
start0=1
for(ii in 1:L)
{
	print(L-ii)
	temp=fread(paste0('G360/',patT[ii],'.txt.gz'),data.table=F)
	#mutemp=fread(paste0('G360_PC/ph12/mu0_',patT[ii],'.txt'),data.table=F)[[1]]##########
	rownames(temp)=temp[,1];temp=temp[rownames(temp1)[INDD],-1]
	pp=dim(temp)[2]
	for(ii in 1:dim(temp)[2])
	{
	med0=median(temp[,ii],na.rm=T)
	mad0=mad(temp[,ii],na.rm=T)
	temp[which(abs(temp[,ii]-med0)>5*mad0),ii]=NA
	indtemp=which(is.na(temp[,ii]))
	if(length(indtemp)>0){#print(length(indtemp));
	temp[indtemp,ii]=median(temp[,ii],na.rm=T)}
	}
	ptemp=dim(temp)[2]
	M1[rownames(temp),start0:(start0+ptemp-1)]=as.matrix(temp)
	#Mu[start0:(start0+ptemp-1)]=mutemp  ##############
	colnames(M1)[start0:(start0+ptemp-1)]=paste0(patT[ii],'_',colnames(temp))
	start0=start0+ptemp
	print(64620-start0)
}
#for(ii in 1:dim(M1)[1]){indtemp=which(is.na(M1[ii,]));M1[ii,indtemp]=Mu[indtemp]} #########
#
M=matrix(NA,dim(temp1)[1],64620)
rownames(M)=rownames(temp1);colnames(M)=1:64620
start0=1
for(ii in 1:L)
{
	print(L-ii)
	temp=fread(subT[ii],data.table=F)
	rownames(temp)=temp[,1];temp=temp[,-1]
	ptemp=dim(temp)[2]
	M[rownames(temp),start0:(start0+ptemp-1)]=as.matrix(temp)
	colnames(M)[start0:(start0+ptemp-1)]=paste0(patT[ii],'_',colnames(temp))
	start0=start0+ptemp
	print(64620-start0)
}

for(ii in 1:dim(M)[1]){if(ii%%1000==0)print(dim(M)[1]-ii);M[ii,]=scale(M[ii,]);}
for(ii in 1:dim(M1)[1]){if(ii%%1000==0)print(dim(M1)[1]-ii);M1[ii,]=scale(M1[ii,]);}
M1=t(M1)
C=M[,1:1000]%*%M1[1:1000,]
NN=ceiling(dim(M)[2]/1000)
for(ii in 2:NN)
{
print(dim(M)[2]-ii*1000)
ind000=((ii-1)*1000+1):min(ii*1000,dim(M)[2])
C=C+M[,ind000]%*%M1[ind000,]
}
rownames(C)=rownames(M)
colnames(C)=colnames(M1)
system(paste0('mkdir -p G360_GlobalEval_out/'))
write.table(C,file=paste0('G360_GlobalEval_out/',iii,'.txt'),quote=F)