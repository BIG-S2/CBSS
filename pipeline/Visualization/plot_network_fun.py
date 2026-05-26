#plot reproducibility along with fiber
#'Auditory','Cingulo-Opercular','Default','Dorsal-Attention','Frontoparietal','Language','Orbito-Affective','Posterior-Multimodal','Somatomotor','Sub','Ventral-Multimodal','Visual1','Visual2'
#cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC;module load python/3.12.1;python3
'''
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC
module load python/3.12.1
python3 plot_network_fun.py Visual1 0 Figure1
python3 plot_network_fun.py Visual2 0 Figure1
python3 plot_network_fun.py Auditory 0 Figure1
python3 plot_network_fun.py Cingulo-Opercular 0 Figure1
python3 plot_network_fun.py Default 0 Figure1
python3 plot_network_fun.py Dorsal-Attention 0 Figure1
python3 plot_network_fun.py Frontoparietal 0 Figure1
python3 plot_network_fun.py Language 0 Figure1
python3 plot_network_fun.py Orbito-Affective 0 Figure1
python3 plot_network_fun.py Posterior-Multimodal 0 Figure1
python3 plot_network_fun.py Somatomotor 0 Figure1
python3 plot_network_fun.py Sub 0 Figure1
python3 plot_network_fun.py Ventral-Multimodal 0 Figure1
python3 plot_network_fun.py Visual1 1 Figure1
python3 plot_network_fun.py Visual2 1 Figure1
python3 plot_network_fun.py Auditory 1 Figure1
python3 plot_network_fun.py Cingulo-Opercular 1 Figure1
python3 plot_network_fun.py Default 1 Figure1
python3 plot_network_fun.py Dorsal-Attention 1 Figure1
python3 plot_network_fun.py Frontoparietal 1 Figure1
python3 plot_network_fun.py Language 1 Figure1
python3 plot_network_fun.py Orbito-Affective 1 Figure1
python3 plot_network_fun.py Posterior-Multimodal 1 Figure1
python3 plot_network_fun.py Somatomotor 1 Figure1
python3 plot_network_fun.py Sub 1 Figure1
python3 plot_network_fun.py Ventral-Multimodal 1 Figure1
'''
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
net0=sys.argv[1];if_within=int(sys.argv[2]);out0=sys.argv[3]
#net0='Visual2';if_within=1;out0='Figure1'
def distinct_hsl(n=100):
    out = []
    phi = (5**0.5 - 1) / 2  # golden-ratio conjugate for spacing
    h = 0.0
    for i in range(n):
        h = (h + phi) % 1.0
        # alternate lightness/saturation to increase separability
        l = 0.55 if (i % 2 == 0) else 0.45
        s = 0.85 if (i % 3) else 0.75
        r, g, b = colorsys.hls_to_rgb(h, l, s)  # note: HLS order
        out.append(f"rgb({int(r*255)},{int(g*255)},{int(b*255)})")
    return out

def add_nifti_isosurface(fig, nii_path, level=0.5, opacity=0.12, color='lightgray', step_size=2):
    img = nb.load(nii_path)
    vol = img.get_fdata()
    mask = (vol > 0)
    verts_vox, faces, _, _ = marching_cubes(mask.astype(float), level=level, step_size=step_size)
    verts_mm = apply_affine(img.affine, verts_vox)
    x, y, z = verts_mm.T
    i, j, k = faces.T
    fig.add_trace(go.Mesh3d(
        x=x, y=y, z=z,
        i=i, j=j, k=k,
        color=color,
        opacity=opacity,
        flatshading=True,
        name='Brain surface',
        showscale=False
    ))

#ROI=sys.argv[1]
tractogram_file = '../template/atlas_100clust_100points.trk'
reference_anatomy_file = '../template/ENIGMA_DTI_FA.nii.gz'
tractogram = load_tractogram(tractogram_file, reference_anatomy_file)
st = StatefulTractogram(tractogram.streamlines, reference_anatomy_file, space=Space.RASMM)
num_points = 100
fa_values_file = '../template/Fiber_atlas_100points_FullName_qc.csv'#'reprod.csv'#'ICC.csv'
#fa_values = np.loadtxt(fa_values_file,delimiter=',')
fa_values = pd.read_csv(fa_values_file, na_values=['NA']) #np.flatnonzero(fa_values['Clust_Dis20mm'].isna())
Utemp=st.streamlines
indd=np.flatnonzero(fa_values['Clust_Dis20mm'].notna())
#indd = np.flatnonzero(fa_values['Clust_Dis20mm'].notna() &(fa_values['Net1_Name'].fillna('')!= fa_values['Net2_Name'].fillna('')))
if if_within==1:
    indd = np.flatnonzero(fa_values['Clust_Dis20mm'].notna() &(fa_values['Net1_Name'].fillna('')==net0) & (fa_values['Net2_Name'].fillna('')==net0))
    #indd11 = np.flatnonzero(fa_values['Clust_Dis20mm'].notna() &(fa_values['Net1_Name'].fillna('')==net0) & (fa_values['Net2_Name'].fillna('')==net0))

if if_within==0:
    indd1 = np.flatnonzero(fa_values['Clust_Dis20mm'].notna() &(fa_values['Net1_Name'].fillna('')==net0) & (fa_values['Net2_Name'].fillna('')!=net0))
    indd2 = np.flatnonzero(fa_values['Clust_Dis20mm'].notna() &(fa_values['Net2_Name'].fillna('')==net0) & (fa_values['Net1_Name'].fillna('')!=net0))
    indd=np.union1d(indd1,indd2)
    #indd11=np.union1d(indd1,indd2)

fa_values, cats0 = pd.factorize(fa_values['Network_Pair'].astype(str).str.strip())
#fa_values1=fa_values.copy();fa_values1[:]=99;fa_values1[indd11]=fa_values[indd11];fa_values=fa_values1.copy()
fa_values=np.repeat(fa_values.reshape(-1, 1), 100, axis=1)
Utemp=Utemp[indd];fa_values=fa_values[indd,:];
streamlines = set_number_of_points(Utemp, num_points)
fig = go.Figure()
ii=0
all_x, all_y, all_z, all_colors, all_opacity = [], [], [], [], []
for idx, streamline in enumerate(streamlines):
    x, y, z = streamline.T
    all_x.extend(x)
    all_y.extend(y)
    all_z.extend(z)
    fa_line = fa_values[idx]
    norm = fa_line
    opacity = np.clip(norm, 0, 1).astype('float')
    opacity[:]=1.0
    if norm[0]==99:
        opacity[:]=0.01
    all_opacity.append(opacity) 
    cmap = distinct_hsl(100)
    cmap[99]="rgb(128,128,128)"
    #colors = ['#000000'] * len(x) #[cmap[value] for value in norm]
    colors = [cmap[value] for value in norm]
    all_colors.extend(colors)

fig.add_trace(go.Scatter3d(
    x=all_x,
    y=all_y,
    z=all_z,
    mode='markers',
    marker=dict(
        size=1,
        color=all_colors,
        opacity=1,
        showscale=False
    )
))
brain_mask_path = '/overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/template/ENIGMA_DTI_FA.nii.gz'
add_nifti_isosurface(fig, brain_mask_path, level=0.5, opacity=0.10, color='lightgray', step_size=2)


fig.update_layout(scene=dict(
    xaxis=dict(visible=False),
    yaxis=dict(visible=False),
    zaxis=dict(visible=False)
))
if if_within==1:
    out1=out0+'/'+net0+'_'+'within.html'

if if_within==0:
    out1=out0+'/'+net0+'_'+'between.html'

os.system('mkdir -p '+out0)
fig.write_html(out1)

cam_coronal  = dict(eye=dict(x=0,  y=2.0, z=0))   # front (A→P)
cam_sagittal = dict(eye=dict(x=-2.0,y=0,  z=0))   # side (L→R)
cam_axial    = dict(eye=dict(x=0,  y=0,  z=2.0)) # top (S→I)
fig.update_layout(scene_camera=cam_coronal)
pio.write_image(fig, f"{out0}/{net0}_coronal.svg",  width=1600, height=1200, scale=2)
fig.update_layout(scene_camera=cam_sagittal)
pio.write_image(fig, f"{out0}/{net0}_sagittal.svg", width=1600, height=1200, scale=2)
fig.update_layout(scene_camera=cam_axial)
pio.write_image(fig, f"{out0}/{net0}_axial.svg",    width=1600, height=1200, scale=2)
#Image.open(f"{out0}/{net0}_axial.pdf").rotate(-90, expand=True).save(f"{out0}/{net0}_axial.pdf")