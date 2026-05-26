#Calculate the fiber length and curvation at each point.
#cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/template;module load python/3.12.1;python3
import numpy as np
import nibabel as nb
from PIL import Image
from dipy.io.streamline import load_tractogram
from dipy.io.stateful_tractogram import StatefulTractogram, Space
from dipy.tracking.streamline import set_number_of_points
import plotly.graph_objs as go
import plotly.express as px
import os
import sys
import pandas as pd
import plotly.express as px
import colorsys
from skimage.measure import marching_cubes
from nibabel.affines import apply_affine
import plotly.io as pio
from scipy.interpolate import splprep, splev
tractogram_file = '../template/atlas_100clust_100points.trk'
reference_anatomy_file = '../template/ENIGMA_DTI_FA.nii.gz'
tractogram = load_tractogram(tractogram_file, reference_anatomy_file)
st = StatefulTractogram(tractogram.streamlines, reference_anatomy_file, space=Space.RASMM)
num_points = 100
fa_values_file = '../template/Fiber_atlas_100points_FullName_qc.csv'#'reprod.csv'#'ICC.csv'
#fa_values = np.loadtxt(fa_values_file,delimiter=',')
fa_values = pd.read_csv(fa_values_file, na_values=['NA']) #np.flatnonzero(fa_values['Clust_Dis20mm'].isna())
Utemp=st.streamlines
Utemp1=np.asarray(Utemp)
FiberLen=np.zeros(Utemp1.shape[0])/0
Curvature=np.zeros(Utemp1.shape[0])/0
t = np.linspace(0, 1, 100)
L=Utemp1.shape[0]
for ii in np.arange(L):
 if np.mod(ii,100)==0:
  print(L-ii)
 coords=Utemp1[ii,:,:]
 diffs = np.diff(coords, axis=0)
 segment_lengths = np.linalg.norm(diffs, axis=1)
 FiberLen[ii] = np.sum(segment_lengths)
 tck, _ = splprep(coords.T, s=0.01)
 r1, r2 = splev(t, tck, der=1), splev(t, tck, der=2)
 r1, r2 = np.array(r1).T, np.array(r2).T
 curvature = np.linalg.norm(np.cross(r1, r2), axis=1) / (np.linalg.norm(r1, axis=1)**3)
 Curvature[ii]=np.mean(curvature)

df = pd.read_csv("Fiber_atlas_100points_FullName_qc.csv")
df["FiberLen"] = FiberLen
df["Curvature"] = Curvature
df.to_csv("Fiber_atlas_100points_FullName_qc_LenCurv.csv", index=False)