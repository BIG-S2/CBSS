'''
#https://stackoverflow.com/questions/71160423/how-to-sample-points-in-3d-in-python-with-origin-and-normal-vector
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/;module purge
export PATH=/proj/tengfei/software/Anaconda/Anaconda3/bin:${PATH}
conda activate
#python 
python3 fibskel_FA_step1.py template/atlas_100clust_100points.trk UKB_Retest/1000240/masked_FA.nii.gz step1_FA_UKB_Phase12/1_fiber 1 template/save_coord/ 2
'''
import sys
import numpy as np
import nibabel as nb
import dipy
from dipy.io.streamline import load_tractogram, save_tractogram 
from dipy.tracking._utils import (_mapping_to_voxel, _to_voxel_coordinates)
from dipy.tracking.streamline import orient_by_streamline
from nibabel.processing import smooth_image
'''
template='/overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality//template/atlas_100clust_100points.trk'
fmap='/overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality//UKB_Phase3/1179894/masked_FA.nii.gz'
out='/overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality//step1_FA_UKB_Phase3/1179894_fib'
DL=1 
coorfile='/overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/template/save_coord'
fwhm=2
'''
template=sys.argv[1] #
fmap=sys.argv[2]
out=sys.argv[3]
DL=sys.argv[4]
coorfile=sys.argv[5]
fwhm=float(sys.argv[6])
A=load_tractogram(template,'same')
A1=A._tractogram._streamlines
ll=len(A1)
NN=100
FAproj=np.zeros([ll,NN])/0
affine=A._affine
lin_T, offset = _mapping_to_voxel(affine)
fa_img = nb.load(fmap)
smoothed_fa_img=smooth_image(fa_img, fwhm)
fa_vol = smoothed_fa_img.get_fdata()
for i in range(ll):
	if np.mod(i,100)==0:
		print(ll-i)
	temp=fa_img.get_fdata().astype('int')
	temp[:]=0;coorT=np.load(coorfile+'/'+str(DL)+'_coord'+str(i)+'.npy');L0=A1[i].shape[0]
	for j in np.arange(NN):
		coor=coorT[j,:,:];indtemp=np.where(~np.isnan(coor[:,0]))[0];coor=coor[indtemp,:].astype('int')
		ii, jj, kk = coor.T/40;
		ii1=np.floor(ii).astype('int');jj1=np.floor(jj).astype('int');kk1=np.floor(kk).astype('int'); #start trilinear interpolation
		dd1=ii-ii1;dd2=jj-jj1;dd3=kk-kk1
		c000=fa_vol[ii1, jj1, kk1];c100=fa_vol[ii1+1, jj1, kk1];c010=fa_vol[ii1, jj1+1, kk1];c001=fa_vol[ii1, jj1, kk1+1]
		c101=fa_vol[ii1+1, jj1, kk1+1];c011=fa_vol[ii1, jj1+1, kk1+1];c110=fa_vol[ii1+1, jj1+1, kk1];c111=fa_vol[ii1+1, jj1+1, kk1+1]
		FAs = c000*(1-dd1)*(1-dd2)*(1-dd3)+c100*dd1*(1-dd2)*(1-dd3)+c010*(1-dd1)*dd2*(1-dd3)+c001*(1-dd1)*(1-dd2)*dd3;
		FAs=FAs+c101*dd1*(1-dd2)*dd3+c011*(1-dd1)*dd2*dd3+c110*dd1*dd2*(1-dd3)+c111*dd1*dd2*dd3;
		FAproj[i,j]=np.max(FAs)

np.save(out+'.npy',FAproj)
