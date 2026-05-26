clear 
home_dir = '/overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/';
%Mat='step2_FA_UKB_Phase5';
Mat='step2_FA_ADNIGO2';
%Mat='step2_FA_ADNI';
code_dir = sprintf('%s/code_%s',home_dir,Mat);
mkdir(code_dir); 
mkdir(sprintf('%s/%s',home_dir,Mat));
system(sprintf('rm %s/*',code_dir))
fid0=fopen(sprintf('%s/All.sh',code_dir),'w'); 
fprintf(fid0,'#!/bin/bash \n');
nn=6090;
K=609;
for iii=1:K
if exist(sprintf('%s/%s/%i_Aligned.mat',home_dir,Mat,iii))
continue;
end
disp(iii)
batNames = sprintf('%s/bat%i.pbs',code_dir,iii);
fid = fopen(batNames,'w');
fprintf(fid,'#!/bin/bash \n');
fprintf(fid,'#SBATCH --ntasks=1 \n');
fprintf(fid,'#SBATCH --time=00:29:59 \n');%1 hour for ukb phase 4; 29 min for hcp
fprintf(fid,'#SBATCH --mem=12000 \n');%12gb for HCP; 24 for ukb phase4
fprintf(fid,'#SBATCH --wrap=mytask \n');
fprintf(fid,'cd %s \n',code_dir);
fprintf(fid,'module load matlab \n');
fprintf(fid,'matlab -nodesktop -nosplash -nodisplay -r ''bat%i'' -logfile bat%i.out \n',iii,iii); 
fclose(fid); 
fprintf(fid0,'sbatch bat%i.pbs\n',iii);
batNames = sprintf('%s/bat%i.m',code_dir,iii);
fid = fopen(batNames,'w');
fprintf(fid,'clear all; \n');
fprintf(fid,'close all; \n');
fprintf(fid,'addpath(genpath(''../elastic'')) \n');
fprintf(fid,'Mat = ''%s''; \n',Mat);
fprintf(fid,'iii=1;\n');
fprintf(fid,'iii=%i \n',iii);
fprintf(fid,'N=10; \n');
fprintf(fid,'lambdaT=table2array(readtable(''%s/plot_and_QC/herit_lambda/HeritOptGWAS.txt''));\n',home_dir);
fprintf(fid,'VoxelT=((iii-1)*N+1):min(iii*N,%i);\n',nn);
fprintf(fid,'N1=%i;iii1=ceil(VoxelT(1)/N1);iii2=VoxelT(1)-(iii1-1)*N1; \n',6090);%%%Change K to 6090
fprintf(fid,'fT=h5read(sprintf(''../%%s/%%i.h5'',Mat,iii1),[''/'' ''M'']);\n');
fprintf(fid,'iii3=min(iii2+N-1,size(fT,2));fT=fT(:,iii2:iii3,:);\n');
fprintf(fid,'if (size(fT,2)<N) & (iii1<10)\n');
fprintf(fid,'fT1=h5read(sprintf(''../%%s/%%i.h5'',Mat,iii1+1),[''/'' ''M'']);fT1=fT1(:,1:(N-size(fT,2)),:);fT=cat(2,fT,fT1);\n');
fprintf(fid,'end\n');
fprintf(fid,'fT=permute(fT,[1,3,2]);\n');
fprintf(fid,'f_aligned = zeros([size(fT,2),size(fT,1),N])/0; \n');
fprintf(fid,'fmeanT=dlmread(sprintf(''../atlas/%%i_Mean.txt'',iii)); \n');
%fprintf(fid,'L21=zeros(size(fT,1),length(VoxelT))/0; \n');
fprintf(fid,'L21=f_aligned; \n');
fprintf(fid,'TT=size(fT,2); \n');
fprintf(fid,'for jjj=1:length(VoxelT)  \n');
fprintf(fid,'	disp(length(VoxelT)-jjj)  \n');
fprintf(fid,'	f=fT(:,:,jjj)'';  \n');
fprintf(fid,'	M1=size(f,1);N11=size(f,2); \n'); %%
fprintf(fid,'	t=(1:TT)/TT;  \n');
fprintf(fid,'	fmean1=fmeanT(:,jjj);  \n');
fprintf(fid,'	a = size(t,1);  \n');
fprintf(fid,'	if (a ~=1)  \n');
fprintf(fid,'		t = t'';  \n');
fprintf(fid,'	end  \n');
fprintf(fid,'	binsize = mean(diff(t));  \n');
fprintf(fid,'	f0 = [fmean1,f];  \n');
fprintf(fid,'	[~,fy] = gradient(f0,binsize,binsize);  \n');
fprintf(fid,'	q = fy./sqrt(abs(fy)+eps);  \n');
fprintf(fid,'	q1 = q(:,1)'';  \n');
fprintf(fid,'	f_1 = f0(:,1);  \n');
fprintf(fid,'	for ii=2:size(f0,2)  \n');
fprintf(fid,'	if mod(ii,100)==0  \n');
fprintf(fid,'	disp(size(f0,2)-ii)  \n');
fprintf(fid,'	end  \n');
fprintf(fid,'		q2 = q(:,ii)'';  \n');
fprintf(fid,'		f_2 = f0(:,ii);  \n');
fprintf(fid,'		[G,T] = DynamicProgrammingQ2_1(q1/norm(q1),t,q2/norm(q2),t,t,t,lambdaT(VoxelT(jjj),1));    \n');
fprintf(fid,'		gam0 = interp1(T,G,t);  \n');
fprintf(fid,'		gam_f = (gam0-gam0(1))/(gam0(end)-gam0(1));  \n');
fprintf(fid,'		f_aligned(:,ii-1,jjj) = interp1(t, f_2, (t(end)-t(1)).*gam_f + t(1))'';  \n');
fprintf(fid,'	    L21(:,ii-1,jjj)=gam_f;  \n');		
fprintf(fid,'	end	  \n');
%fprintf(fid,'	L21(:,jjj)=sum((f_aligned(:,:,jjj)-fmean1*ones(1,size(f_aligned,2))).^2);  \n');
fprintf(fid,'	%%med1=median(L21(:,jjj));mad1=mad(L21(:,jjj));indqc1=(abs(L21(:,jjj)-med1)>5*mad1);  \n');
fprintf(fid,'end  \n');
fprintf(fid,'save(sprintf(''../%%s/%%i_%%s.mat'',Mat,iii,''Aligned''),''f_aligned'');  \n');
fprintf(fid,'save(sprintf(''../%%s/%%i_%%s.mat'',Mat,iii,''TimeShift''),''L21'');  \n');
fprintf(fid,'clear;exit;return;  \n');
fclose(fid); 
end; 
fclose(fid0); 
clear all; 
