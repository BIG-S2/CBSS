#This script will randomly pick 300 subjects for the use of generating the Karcher Mean.
#cd step2_FA_UKB_Retest/ #change this path to your own output folder
#module load r/4.1.0;R
library(rhdf5)
library(R.matlab)
A=h5read('1.h5','M')
A=A[1:300,,]
system(paste0('mkdir ../atlas/'))
writeMat('/overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/atlas/atlas.mat',fT=aperm(A,c(1,3,2)))
