cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC/fMRI_coupling/PopulationAnalysis/
module load r/4.4.0;R
library(data.table)
Input='G360_CBSS_Line_mask_out'#'G360_out'
Out=gsub('G360','G360_GlobalEval',Input)
subT=Sys.glob(paste0(Input,'/*_pred.txt'))
patT=gsub(paste0(Input,'/'),'',gsub('_pred.txt','',subT))
L=length(patT)
temp1=fread(subT[1],data.table=F)
rownames(temp1)=temp1[,1];temp1=temp1[,-1]
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
##read in the true data
M1=matrix(NA,dim(temp1)[1],64620)
rownames(M1)=rownames(temp1);colnames(M1)=1:64620
start0=1
for(ii in 1:L)
{
	print(L-ii)
	temp=fread(paste0('G360/',patT[ii],'.txt.gz'),data.table=F)
	rownames(temp)=temp[,1];temp=temp[rownames(temp1),-1]
	ptemp=dim(temp)[2]
	M1[rownames(temp),start0:(start0+ptemp-1)]=as.matrix(temp)
	colnames(M1)[start0:(start0+ptemp-1)]=paste0(patT[ii],'_',colnames(temp))
	start0=start0+ptemp
	print(64620-start0)
}
mu_tr=rep(NA,64620)
start0=1
for(ii in 1:L)
{
	print(L-ii)
	temp=read.table(paste0('G360_PC/ph12/mu0_',patT[ii],'.txt'))[[1]]
	ptemp=length(temp)
	mu_tr[start0:(start0+ptemp-1)]=as.numeric(temp)
	start0=start0+ptemp
}
for(jj in 1:dim(M)[2]){if(jj%%1000==0)print(dim(M)[2]-jj);M[,jj]=M[,jj]-mu_tr[jj];M1[,jj]=M1[,jj]-mu_tr[jj];}
S=rep(NA,dim(M)[1])
for(ii in 1:dim(M)[1]){
	if(ii%%200==0){print(dim(M)[1]-ii);print(mean(S[1:(ii-1)]))}
	S[ii]=cor.test(M[ii,],M1[ii,])$est
}
names(S)=rownames(M)
write.table(S,file=paste0(Out,'/globalCentralCor.txt'),quote=F)
#prediction correlation: 0.146 # for masked CBSS_line:  
#############################################################
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC/fMRI_coupling/PopulationAnalysis/
module load r/4.4.0;R
library(data.table)
subT=Sys.glob(paste0('G360_out/*_pred.txt'))
patT=gsub('G360_out/','',gsub('_pred.txt','',subT))
L=length(patT)
temp1=fread(subT[1],data.table=F)
rownames(temp1)=temp1[,1];temp1=temp1[,-1]
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
##read in the true data
M1=matrix(NA,dim(temp1)[1],64620)
rownames(M1)=rownames(temp1);colnames(M1)=1:64620
start0=1
for(ii in 1:L)
{
	print(L-ii)
	temp=fread(paste0('G360/',patT[ii],'.txt.gz'),data.table=F)
	rownames(temp)=temp[,1];temp=temp[rownames(temp1),-1]
	ptemp=dim(temp)[2]
	M1[rownames(temp),start0:(start0+ptemp-1)]=as.matrix(temp)
	colnames(M1)[start0:(start0+ptemp-1)]=paste0(patT[ii],'_',colnames(temp))
	start0=start0+ptemp
	print(64620-start0)
}
for(ii in 1:dim(M)[1]){if(ii%%1000==0)print(dim(M)[1]-ii);M[ii,]=scale(M[ii,]);M1[ii,]=scale(M1[ii,]);}
C=matrix(NA,dim(M1)[1],dim(M)[1])
M1=t(M1)
for(ii in 1:dim(C)[1])
{
print(dim(M)[1]-ii)
C[1:4000,]=M[1:1000,1:1000]%*%M1[1:1000,]
}

cb_metrics <- function(M, M1) {
  M <- as.matrix(M) #A_pred; A_pred1;A_test1
  M1 <- as.matrix(M1)
  n <- nrow(M1)
  row_corr <- function(A, B) {
    A <- scale(A, center = TRUE, scale = FALSE)
    B <- scale(B, center = TRUE, scale = FALSE)
    num <- rowSums(A * B)
    den <- sqrt(rowSums(A^2) * rowSums(B^2))
    num / den
  }
  r_rows_demean <- row_corr(sweep(M, 2, mu_tr),sweep(M1, 2, mu_tr))
  avgcorr_demean <- mean(r_rows_demean, na.rm = TRUE)#0.157
  #3. avgrank
  Xn <- t(scale(t(M1)))
  Yn <- t(scale(t(M)))
  C <- Xn %*% t(Yn)/dim(Xn)[2]  # n×n correlation matrix
  mm=rep(NA,dim(Xn)[1]);for(ii in 1:dim(Xn)[1])mm[ii]= mean(C[ii,]<C[ii,ii])
  avgrank <- mean(mm) #0.654
  avgcorr <- mean(diag(C))
  data.frame(avgcorr, avgcorr_demean, avgrank)
}
cb_metrics0=cb_metrics(M, M1, A_train1)
write.table(cb_metrics0,file=paste0('G360_out/',outfile,'_metrics.txt'),quote=F)
write.table(M,file=paste0('G360_out/',outfile,'_pred.txt'),quote=F)
#############################################################
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC/fMRI_coupling/PopulationAnalysis/
module load r/4.4.0;R
library(data.table)
Input0='G360_CBSS_Line_mask_out'#'G360_out'
Input=gsub('G360','G360_GlobalEval',Input0)
C=matrix(NA,18095,18095)
for(ii in 1:61)
{
	print(61-ii)
	A0=fread(paste0(Input,'/',ii,'.txt'),data.table=F)
	rownames(A0)=A0[,1];A0=A0[,-1]
	indd=((ii-1)*300+1):min((ii*300),18095)
	C[indd,]=t(as.matrix(A0))
}
indd0=which(!is.na(C[,1]))
C=C[indd0,indd0]
S=rep(NA,dim(C)[1]);for(ii in 1:dim(C)[1])S[ii]=mean(C[ii,-ii]<C[ii,ii])
avgrank <- mean(S) #0.896
avgcorr <- mean(diag(C))

mu_tr <- colMeans(A_train1)
r_rows_demean <- row_corr(sweep(A_pred1, 2, mu_tr),sweep(A_test1, 2, mu_tr))
avgcorr_demean <- mean(r_rows_demean, na.rm = TRUE)#0.157