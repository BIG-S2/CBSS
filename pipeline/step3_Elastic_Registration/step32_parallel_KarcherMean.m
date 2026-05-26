%This file will submit parallel jobs to generate Karcher mean for all 6090 pathways among the population.
clear 
home_dir = '/overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/' 
code_dir = fullfile(home_dir,'code1'); 
system(sprintf('mkdir -p %s/atlas/',home_dir))
mkdir(code_dir); 
delete(sprintf('%s/*',code_dir)); 
fid0=fopen(sprintf('%s/All.sh',code_dir),'w'); 
fprintf(fid0,'#!/bin/bash \n');
nn=6090;
K=609;
for iii=1:K  %
if exist(sprintf('%s/atlas/%i_Aligned.mat',home_dir,iii))
continue;
end
batNames = sprintf('%s/bat%i.pbs',code_dir,iii);
fid = fopen(batNames,'w');
fprintf(fid,'#!/bin/bash \n');
fprintf(fid,'#SBATCH --ntasks=1 \n');
fprintf(fid,'#SBATCH --time=3:59:59 \n');%10 for 4:00:00
fprintf(fid,'#SBATCH --mem=12000 \n');
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
fprintf(fid,'lambda = 0; \n');
fprintf(fid,'load ../atlas/atlas.mat \n');
fprintf(fid,'TT=size(fT,2);\n');
fprintf(fid,'iii=1;%%1:610 \n');
fprintf(fid,'iii=%i \n',iii);
fprintf(fid,'N=10; \n');
fprintf(fid,'VoxelT=((iii-1)*N+1):min(iii*N,%i);%%73 \n',nn);
fprintf(fid,'f_aligned = zeros([size(fT,2),size(fT,1),N])/0; \n');
fprintf(fid,'fmeanT=zeros(size(fT,2),length(VoxelT))/0; \n');
fprintf(fid,'L21=zeros(size(fT,1),length(VoxelT))/0; \n');
fprintf(fid,'for jjj=1:length(VoxelT) \n');
fprintf(fid,'	disp(length(VoxelT)-jjj) \n');
fprintf(fid,'	f=fT(:,:,VoxelT(jjj))''; \n');
fprintf(fid,'	M1=size(f,1);N1=size(f,2); \n');
fprintf(fid,'	for r = 1:4 \n');
fprintf(fid,'		for i = 1:N1 \n');
fprintf(fid,'			f(2:(M1-1),i) = (f(1:(M1-2),i)+2*f(2:(M1-1),i) + f(3:M1,i))/4; \n');
fprintf(fid,'		end \n');
fprintf(fid,'	end \n');
fprintf(fid,'	t=(1:TT)/TT; \n');
fprintf(fid,'	[fn,~,~,fmean]=time_warping(f,t,lambda); \n');
fprintf(fid,'	L2=sum((fn-fmean*ones(1,size(fn,2))).^2); \n');
fprintf(fid,'	med0=median(L2);mad0=mad(L2); \n');
fprintf(fid,'	indqc=find(abs(L2-med0)>5*mad0); \n');
fprintf(fid,'	[fn1,~,~,fmean1]=time_warping(f(:,setdiff(1:size(f,2),indqc)),t,lambda); \n');
fprintf(fid,'	fmeanT(:,jjj)=fmean1; \n');
fprintf(fid,'	a = size(t,1); \n');
fprintf(fid,'	if (a ~=1) \n');
fprintf(fid,'		t = t''; \n');
fprintf(fid,'	end \n');
fprintf(fid,'	binsize = mean(diff(t)); \n');
fprintf(fid,'	f0 = [fmean1,f]; \n');
fprintf(fid,'	[~,fy] = gradient(f0,binsize,binsize); \n');
fprintf(fid,'	q = fy./sqrt(abs(fy)+eps); \n');
fprintf(fid,'	q1 = q(:,1)''; \n');
fprintf(fid,'	f_1 = f0(:,1); \n');
fprintf(fid,'	for ii=2:size(f0,2) \n');
fprintf(fid,'		q2 = q(:,ii)''; \n');
fprintf(fid,'		f_2 = f0(:,ii); \n');
fprintf(fid,'		[G,T] = DynamicProgrammingQ2(q1/norm(q1),t,q2/norm(q2),t,t,t); \n');
fprintf(fid,'		gam0 = interp1(T,G,t); \n');
fprintf(fid,'		gam_f = (gam0-gam0(1))/(gam0(end)-gam0(1)); \n');
fprintf(fid,'		f_aligned(:,ii-1,jjj) = interp1(t, f_2, (t(end)-t(1)).*gam_f + t(1))''; \n');
fprintf(fid,'	end	 \n');
fprintf(fid,'	L21(:,jjj)=sum((f_aligned(:,:,jjj)-fmean1*ones(1,size(f_aligned,2))).^2); \n');
fprintf(fid,'	%%med1=median(L21);mad1=mad(L21);indqc1=(abs(L21-med1)>5*mad1); \n');
fprintf(fid,'end \n');
fprintf(fid,'save(sprintf(''../atlas/%%i_%%s.mat'',iii,''Aligned''),''f_aligned''); \n');
fprintf(fid,'dlmwrite(sprintf(''../atlas/%%i_%%s.txt'',iii,''Mean''),fmeanT); \n');
fprintf(fid,'dlmwrite(sprintf(''../atlas/%%i_%%s.txt'',iii,''L2''),L21); \n');
fclose(fid); 
end; 
fclose(fid0); 
clear all; 
