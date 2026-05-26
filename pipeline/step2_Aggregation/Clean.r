#######################################################################################################################################
#This file will save population-level FA profiles into 10 .h5 files for UKB phase 12, and UKB phase 3 data;
#This file will save population-level FA profiles into one single .h5 files for UKB other phase data and other lifespan datasets.
#######################################################################################################################################
#srun -t 8:00:00 -p interact -N 1 -n 1 --x11=first --mem=16gb --pty /bin/bash
#cd projects/UKB_PSC_CommonAtlas/OtherModality
#module load r/4.1.0
#module load gcc
#R
dir0='step1_FA_UKB_Phase12/'
out='step2_FA_UKB_Phase12/'
system(paste0('mkdir ',out))
sub0=dir(dir0)
write.table(gsub('_fib.npy','',sub0),file=paste0(out,'UKB_Phase12.txt'),quote=F,row.names=F,col.names=F)
L=length(sub0)
library(reticulate)
library(rhdf5)
np=import('numpy')
temp=np$load(paste0(dir0,sub0[1]))
p=dim(temp)[1];N=dim(temp)[2]
NN=ceiling(p/10)
for(ii in 1:10)
{
	Vec=((ii-1)*NN+1):min(ii*NN,p)
	M=array(NA,c(L,length(Vec),N))
	for(i in 1:L)
	{
		if(i%%100==0)print(L-i)
		temp=np$load(paste0(dir0,sub0[i]))
		M[i,,]=round(temp[Vec,],3)
	}
	h5createFile(paste0(out,ii,'.h5'))
	h5write(M, paste0(out,ii,'.h5'), "M")
}
#####################################################################################################################
#srun -t 8:00:00 -p interact -N 1 -n 1 --x11=first --mem=16gb --pty /bin/bash
#cd projects/UKB_PSC_CommonAtlas/OtherModality
#module load r/4.1.0
#module load gcc
#R
dir0='step1_FA_UKB_Phase3/'
out='step2_FA_UKB_Phase3/'
system(paste0('mkdir ',out))
sub0=dir(dir0)
write.table(gsub('_fib.npy','',sub0),file=paste0(out,'UKB_Phase3.txt'),quote=F,row.names=F,col.names=F)
L=length(sub0)
library(reticulate)
library(rhdf5)
np=import('numpy')
temp=np$load(paste0(dir0,sub0[1]))
p=dim(temp)[1];N=dim(temp)[2]
NN=ceiling(p/10)
for(ii in 1:10)
{
	Vec=((ii-1)*NN+1):min(ii*NN,p)
	M=array(NA,c(L,length(Vec),N))
	for(i in 1:L)
	{
	if(i%%100==0)print(L-i)
	temp=np$load(paste0(dir0,sub0[i]))
	M[i,,]=round(temp[Vec,],3)
	}
	#M=round(M,3)
	h5createFile(paste0(out,ii,'.h5'))
	h5write(M, paste0(out,ii,'.h5'), "M")
}
#####################################################################################################################
#srun -t 8:00:00 -p interact -N 1 -n 1 --x11=first --mem=16gb --pty /bin/bash
#cd projects/UKB_PSC_CommonAtlas/OtherModality
#module load r/4.1.0
#module load gcc
#R
ii=1
dir0='step1_FA_UKB_Phase5/';out='step2_FA_UKB_Phase5/' #for Phase4 and Retest data, change "Phase5" to "Phase4" and "Retest"
system(paste0('mkdir ',out))
sub0=dir(dir0)
write.table(gsub('_fib.npy','',sub0),file=paste0(out,'UKB_Phase5.txt'),quote=F,row.names=F,col.names=F)
L=length(sub0)
library(reticulate)
library(rhdf5)
np=import('numpy')
temp=np$load(paste0(dir0,sub0[1]))
p=dim(temp)[1];N=dim(temp)[2]
NN=p
Vec=1:p
M=array(NA,c(L,length(Vec),N))
for(i in 1:L)
{
if(i%%100==0)print(L-i)
temp=np$load(paste0(dir0,sub0[i]))
M[i,,]=round(temp[Vec,],3)
}
h5createFile(paste0(out,ii,'.h5'))
h5write(M, paste0(out,ii,'.h5'), "M")
#####################################################################################################################
#srun -t 8:00:00 -p interact -N 1 -n 1 --x11=first --mem=16gb --pty /bin/bash
#cd projects/UKB_PSC_CommonAtlas/OtherModality
#module load r/4.1.0
#module load gcc
#R
ii=1
dir0='step1_FA_UKB_Retest/'
out='step2_FA_UKB_Retest/'
system(paste0('mkdir ',out))
sub0=dir(dir0)
write.table(gsub('_fib.npy','',sub0),file=paste0(out,'UKB_Retest.txt'),quote=F,row.names=F,col.names=F)
L=length(sub0)
library(reticulate)
library(rhdf5)
np=import('numpy')
temp=np$load(paste0(dir0,sub0[1]))
p=dim(temp)[1];N=dim(temp)[2]
NN=p
Vec=1:p
M=array(NA,c(L,length(Vec),N))
for(i in 1:L)
{
	if(i%%100==0)print(L-i)
	temp=np$load(paste0(dir0,sub0[i]))
	M[i,,]=round(temp[Vec,],3)       
}
h5createFile(paste0(out,ii,'.h5'))
h5write(M, paste0(out,ii,'.h5'), "M")
####################################################################################################################
#srun -t 8:00:00 -p interact -N 1 -n 1 --x11=first --mem=64gb --pty /bin/bash
#cd projects/UKB_PSC_CommonAtlas/OtherModality
#module load r/4.1.0
#module load gcc
#R
ii=1
dataname='PNC'#'HCPA' #HCP
dir0=paste0('step1_FA_',dataname,'/')
out=paste0('step2_FA_',dataname,'/')
system(paste0('mkdir ',out))
sub0=dir(dir0)
write.table(gsub('_fib.npy','',sub0),file=paste0(out,dataname,'.txt'),quote=F,row.names=F,col.names=F)
L=length(sub0)
library(reticulate)
library(rhdf5)
np=import('numpy')
temp=np$load(paste0(dir0,sub0[1]))
p=dim(temp)[1];N=dim(temp)[2]
NN=p
Vec=1:p
M=array(NA,c(L,length(Vec),N))
for(i in 1:L)
{
	if(i%%100==0)print(L-i)
	temp=np$load(paste0(dir0,sub0[i]))
	M[i,,]=round(temp[Vec,],3)
}
h5createFile(paste0(out,ii,'.h5'))
h5write(M, paste0(out,ii,'.h5'), "M")
####################################################################################################################
#srun -t 8:00:00 -p interact -N 1 -n 1 --x11=first --mem=64gb --pty /bin/bash
#cd projects/UKB_PSC_CommonAtlas/OtherModality
#module load r/4.1.0
#module load gcc
#R
ii=1 #1:1
dataname='ADNIGO2'#'PING'#'PNC'
dir0=paste0('step1_FA_',dataname,'/')
out=paste0('step2_FA_',dataname,'/')
system(paste0('mkdir ',out))
sub0=dir(dir0)
write.table(gsub('_fib.npy','',sub0),file=paste0(out,dataname,'.txt'),quote=F,row.names=F,col.names=F)
L=length(sub0)
library(reticulate)
library(rhdf5)
np=import('numpy')
temp=np$load(paste0(dir0,sub0[1]))
p=dim(temp)[1];N=dim(temp)[2]
NN=p
Vec=1:p
M=array(NA,c(L,length(Vec),N))
for(i in 1:L)
{
	if(i%%100==0)print(L-i)
	temp=np$load(paste0(dir0,sub0[i]))
	M[i,,]=round(temp[Vec,],3)
}
h5createFile(paste0(out,ii,'.h5'))
h5write(M, paste0(out,ii,'.h5'), "M")

