'''
#https://stackoverflow.com/questions/71160423/how-to-sample-points-in-3d-in-python-with-origin-and-normal-vector
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/template;module purge
export PATH=/proj/tengfei/software/Anaconda/Anaconda3/bin:${PATH}
conda activate
#python 
python3 fibskel_coord.py atlas_100clust_100points.trk 1
'''
import sys
import numpy as np
import nibabel as nb
import dipy
from dipy.io.streamline import load_tractogram, save_tractogram 
from dipy.tracking._utils import (_mapping_to_voxel, _to_voxel_coordinates)
from dipy.tracking.streamline import orient_by_streamline
import os
os.system('mkdir save_coord/')
'''
template='atlas_100clust_100points.trk';
DL=1
fmap='ENIGMA_DTI_FA.nii.gz'
'''
template=sys.argv[1]
#fmap=sys.argv[2]
#DL=np.float(sys.argv[3])
DL=np.float(sys.argv[2])
A=load_tractogram(template,'same')
A1=A._tractogram._streamlines
ll=len(A1)
NN=60;#or 200
FAproj=np.zeros([ll,NN])/0
affine=A._affine
lin_T, offset = _mapping_to_voxel(affine)
#fa_img = nb.load(fmap)
#fa_vol = fa_img.get_data()
for i in range(ll):
	if np.mod(i,10)==0:
		print(ll-i)
	dif0=np.diff(A1[i],axis=0);dif0=np.concatenate((dif0,dif0[-1,:][np.newaxis,:]));L0=A1[i].shape[0];coorT=np.zeros([L0,600,3])/0
	for j in np.arange(L0):
		vtemp=np.array([0,dif0[j,2],-dif0[j,1]]);utemp=np.cross(vtemp,dif0[j,:])
		temp0=np.arange(NN).reshape([NN,1])/(NN+0.);temp1=vtemp*np.cos(np.pi*2*temp0)+utemp*np.sin(np.pi*2*temp0)
		norm0=np.sqrt(np.sum(temp1**2,axis=1));temp11=1/norm0.reshape([len(norm0),1])*temp1
		M=np.zeros([NN*10,3])/0
		for k in np.arange(10):
			M[(k*NN):((k+1)*NN),:]=A1[i][j,:]+temp11*(DL+0.)/10*k
		M1=_to_voxel_coordinates(M*40, lin_T, offset*40);M1_str=np.asarray(['_'.join(idx for idx in sub) for sub in M1.astype('str')])
		[aa1,ab1]=np.unique(M1_str,return_index=True);coorT[j,:len(ab1),:]=M1[ab1,:]
	np.save('save_coord/'+str(DL)+'_coord'+str(i)+'.npy',coorT)

'''
temp=fa_img.get_data().astype('int')
temp[:]=0
for i in range(ll):
	if np.mod(i,100)==0:
		print(ll-i)
	coorT=np.load('save_coord/'+str(DL)+'_coord'+str(i)+'.npy')
	for j in np.arange(100):
		temp0=coorT[j,:,:];indtemp=np.where(~np.isnan(temp0[:,0]))[0];temp0=temp0[indtemp,:].astype('int')
		[ii,jj,kk]=(temp0.T/40).astype('int')
		temp[ii,jj,kk]=1

img=nb.Nifti1Image(temp,affine)
nb.save(img,'ProjAtlas_'+str(DL)+'.nii.gz')
'''

