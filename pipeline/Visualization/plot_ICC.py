#plot reproducibility along with fiber
#cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot;module load python/3.7.9;python3
import numpy as np
from dipy.io.streamline import load_tractogram
from dipy.io.stateful_tractogram import StatefulTractogram, Space
from dipy.tracking.streamline import set_number_of_points
import plotly.graph_objs as go
import plotly.express as px
import os
import sys
#ROI=sys.argv[1]
tractogram_file = '../template/atlas_100clust_100points.trk'
reference_anatomy_file = '../template/ENIGMA_DTI_FA.nii.gz'
tractogram = load_tractogram(tractogram_file, reference_anatomy_file)
st = StatefulTractogram(tractogram.streamlines, reference_anatomy_file, space=Space.RASMM)
num_points = 100
fa_values_file = 'ICC.csv'#'reprod.csv'#'ICC.csv'
fa_values = np.loadtxt(fa_values_file,delimiter=',')
fa_values_file1 = 'reprod.csv'
fa_values1 = np.loadtxt(fa_values_file1,delimiter=',')
FROI=np.loadtxt('../template/atlas_100clust_100points.csv',delimiter=',',dtype='str')
FROI=FROI[1:,:];ROI1=FROI[:,0].astype('int');ROI2=FROI[:,1].astype('int')
Utemp=st.streamlines
#indd=np.where((np.mean(fa_values,axis=1)>0.65)&(np.mean(fa_values1,axis=1)>0.75))[0]
indd=np.where((np.mean(fa_values,axis=1)>0.65)&(np.mean(fa_values1,axis=1)>0.75)&((ROI1>=347)|(ROI2>=347)))[0]
#indd=np.where(((ROI1>=347)|(ROI2>=347)))[0]
Utemp=Utemp[indd];fa_values=fa_values[indd,:];fa_values1=fa_values1[indd,:]
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
    norm = (fa_line - 0.) / (1 - 0.)
    #norm = (fa_line + 0.25) / (1 + 0.25)
    opacity = np.clip(norm, 0, 1)
    all_opacity.extend(opacity)
    cmap = px.colors.sequential.Rainbow
    colors = [cmap[int(value * (len(cmap) - 1))] for value in norm]
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
        showscale=True,
        colorbar=dict(title='Reliability', titleside='right')
    )
))
fig.update_layout(scene=dict(
    xaxis=dict(showbackground=False, showline=True, linecolor='rgba(0,0,0,0)'),
    yaxis=dict(showbackground=False, showline=True, linecolor='rgba(0,0,0,0)'),
    zaxis=dict(showbackground=False, showline=True, linecolor='rgba(0,0,0,0)')
))
#fig.write_html('reprod.html')
#fig.write_html('ICC.html')
#fig.write_html('reprod_qc.html')
#fig.write_html('ICC_qc.html')
fig.write_html('ICC_subcort_qc.html')
#fig.write_html('ICC_subcort.html')




#############for python/3.12.1
#plot reproducibility along with fiber
#cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/plot_and_QC;module load python/3.12.1;python3
import numpy as np
from dipy.io.streamline import load_tractogram
from dipy.io.stateful_tractogram import StatefulTractogram, Space
from dipy.tracking.streamline import set_number_of_points
import plotly.graph_objs as go
import plotly.express as px
import os
import sys
#ROI=sys.argv[1]
tractogram_file = '../template/atlas_100clust_100points.trk'
reference_anatomy_file = '../template/ENIGMA_DTI_FA.nii.gz'
tractogram = load_tractogram(tractogram_file, reference_anatomy_file)
st = StatefulTractogram(tractogram.streamlines, reference_anatomy_file, space=Space.RASMM)
num_points = 100
fa_values_file = 'ICC.csv'#'reprod.csv'#'ICC.csv'
fa_values = np.loadtxt(fa_values_file,delimiter=',')
fa_values_file1 = 'reprod.csv'
fa_values1 = np.loadtxt(fa_values_file1,delimiter=',')
FROI=np.loadtxt('../template/atlas_100clust_100points.csv',delimiter=',',dtype='str')
FROI=FROI[1:,:];ROI1=FROI[:,0].astype('int');ROI2=FROI[:,1].astype('int')
Utemp=st.streamlines
#indd=np.where((np.mean(fa_values,axis=1)>0.65)&(np.mean(fa_values1,axis=1)>0.75))[0]
indd=np.where((np.mean(fa_values,axis=1)>0.65)&(np.mean(fa_values1,axis=1)>0.75)&((ROI1>=347)|(ROI2>=347)))[0]
#indd=np.where(((ROI1>=347)|(ROI2>=347)))[0]
Utemp=Utemp[indd];fa_values=fa_values[indd,:];fa_values1=fa_values1[indd,:]
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
    norm = (fa_line - 0.) / (1 - 0.)
    #norm = (fa_line + 0.25) / (1 + 0.25)
    opacity = np.clip(norm, 0, 1)
    all_opacity.extend(opacity)
    cmap = px.colors.sequential.Rainbow
    colors = [cmap[int(value * (len(cmap) - 1))] for value in norm]
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
        showscale=True,
        colorbar=dict(title='Reliability', titleside='right')
    )
))


fig.add_trace(go.Scatter3d(
    x=all_x,
    y=all_y,
    z=all_z,
    mode='markers',
    marker=dict(
        size=1,
        color=all_colors,   # numeric array => color scale shown
        opacity=1,
        showscale=True,
        colorbar=dict(title=dict(text='Reliability', side='right'))
    )
))

fig.update_layout(scene=dict(
    xaxis=dict(showbackground=False, showline=True, linecolor='rgba(0,0,0,0)'),
    yaxis=dict(showbackground=False, showline=True, linecolor='rgba(0,0,0,0)'),
    zaxis=dict(showbackground=False, showline=True, linecolor='rgba(0,0,0,0)')
))
#fig.write_html('reprod.html')
#fig.write_html('ICC.html')
#fig.write_html('reprod_qc.html')
#fig.write_html('ICC_qc.html')
fig.write_html('ICC_subcort_qc.html')
#fig.write_html('ICC_subcort.html')

