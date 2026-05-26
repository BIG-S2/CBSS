'''
cd /overflow/tengfei/user/tengfei/projects/UKB_PSC_CommonAtlas/OtherModality/template
module load python/3.12.1;python3
'''
import sys
import numpy as np
import nibabel as nib
from scipy.ndimage import distance_transform_edt


def dilate_labels_by_mm(input_nii, output_nii, dilation_mm=1.0):
    img = nib.load(input_nii)
    data = img.get_fdata()
    affine = img.affine
    header = img.header.copy()
    # Preserve integer labels
    data = np.asarray(data)
    labels = data.astype(np.int32)
    # Get voxel sizes in mm
    voxel_sizes = header.get_zooms()[:3]
    # Original nonzero mask
    orig_mask = labels > 0
    if not np.any(orig_mask):
        raise ValueError("Input image has no nonzero labels.")
    # Distance from every zero voxel to nearest nonzero voxel
    # return_indices gives the coordinates of the nearest source voxel
    distances, indices = distance_transform_edt(
        ~orig_mask,
        sampling=voxel_sizes,
        return_indices=True
    )
    # Define dilated mask: within dilation_mm of original labeled region
    dilated_mask = distances <= dilation_mm
    # Find newly added voxels only
    new_voxels = dilated_mask & (~orig_mask)
    # Initialize output with original labels
    out = labels.copy()
    # Assign each new voxel the label of its nearest original nonzero voxel
    nearest_x = indices[0][new_voxels]
    nearest_y = indices[1][new_voxels]
    nearest_z = indices[2][new_voxels]
    out[new_voxels] = labels[nearest_x, nearest_y, nearest_z]
    # Save result
    out_img = nib.Nifti1Image(out.astype(labels.dtype), affine, header)
    nib.save(out_img, output_nii)
    print(f"Saved dilated label image to: {output_nii}")

dilate_labels_by_mm('JHU-WhiteMatter-labels-1mm.nii.gz', 'JHU-WhiteMatter-labels-1mm_dil.nii.gz', 1.0)
dilate_labels_by_mm('JHU-WhiteMatter-labels-1mm.nii.gz', 'JHU-WhiteMatter-labels-1mm_dil2.nii.gz', 2.0)
dilate_labels_by_mm('JHU-WhiteMatter-labels-1mm.nii.gz', 'JHU-WhiteMatter-labels-1mm_dil3.nii.gz', 3.0)