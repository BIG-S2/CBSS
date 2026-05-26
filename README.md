# **Custom Code for the Image Data Processing and Analysis Pipeline of Connectome-Based Spatial Statistics**

## 1. pipeline/

The pipeline/ folder contains the main processing scripts for constructing fiber-wise and point-wise diffusion representations from subject-level diffusion MRI maps.

pipeline/
├── step0_Neighborhood Coord_1mm/
├── step1_FA_projection/
├── step2_Aggregation/
├── step3_Elastic_Registration/
└── step4_representation/
step0_Neighborhood Coord_1mm/

This step generates local neighborhood coordinates around atlas streamline points.

For each point along each atlas streamline, the script samples a local cross-sectional neighborhood in the plane approximately perpendicular to the streamline direction. The sampled coordinates are converted to voxel space, duplicate neighboring samples are removed, and the resulting coordinate arrays are saved for downstream projection.

These precomputed neighborhood coordinates are later used to extract local maximum diffusion values from subject-specific maps.

step1_FA_projection/

This step projects subject-level diffusion scalar maps, such as FA, MD, RD, L1 or NODDI-derived measures, onto the atlas streamlines.

For each streamline point, the script loads the precomputed neighborhood coordinates, samples the subject-specific diffusion map using trilinear interpolation, and assigns the maximum value within the local neighborhood to the corresponding atlas point.

The output is a fiber-by-point matrix for each subject.

step2_Aggregation/

This step aggregates projected point-wise diffusion values into higher-level representations.

Depending on the analysis goal, this may include summarizing streamline-level values, network-pair values, or region-specific features from the projected fiber profiles.

step3_Elastic_Registration/

This step performs elastic registration of fiber profiles.

Subject-specific fiber profiles are aligned to atlas mean profiles using SRVF-based dynamic programming. The alignment estimates a monotone warping function for each profile, allowing local shifts along the fiber trajectory to be corrected before downstream feature extraction.

The output includes aligned fiber profiles and time-shift/warping functions.

step4_representation/

This step derives final CBSS representations from aligned profiles.

These representations may include profile means, maximum values, PCA scores, or other fiber-wise/network-wise features used for statistical modelling, prediction, heritability analysis, and visualization.

## 2. template/

The template/ folder contains atlas and reference files required by the CBSS pipeline.

template/
├── atlas_100clust_100points_QC.trk
├── atlas_100clust_100points_Total.trk
├── coordinate_atlas.npy
├── Fiber_atlas_100points_FullName_qc_LenCur.csv
├── FiberLen.py
└── roi_id_clean_merged.mat
atlas_100clust_100points_QC.trk

Quality-controlled streamline atlas.

This file contains the final atlas streamlines used for projection and downstream analysis after quality control.

atlas_100clust_100points_Total.trk

Complete streamline atlas before or outside final quality-control filtering.

This file may be used for reference, comparison, or intermediate processing.

coordinate_atlas.npy

NumPy file containing atlas coordinate information.

This file supports coordinate-based indexing, mapping, or projection steps in the pipeline.

Fiber_atlas_100points_FullName_qc_LenCur.csv

Metadata table for atlas fibers.

This file contains fiber names, quality-control information, length/curvature-related information, and possibly network-pair labels used to organize fiber-wise results.

FiberLen.py

Python script for computing or summarizing fiber length-related metrics.

roi_id_clean_merged.mat

MATLAB file containing cleaned and merged ROI identifiers.

This file supports mapping fibers or endpoints to anatomical/functional regions.

## 3. Analysis/

The Analysis/ folder contains downstream statistical analysis and visualization scripts.

Analysis/
├── Aging/
├── Lifespan/
├── Prediction/
├── SC_FC/
└── Visualization/
Aging/

Scripts for aging-related analyses, including age-associated white-matter changes, older-adult effects, and disease-relevant structural-connectivity patterns.

Lifespan/

Scripts for modelling lifespan trajectories of CBSS-derived features.

These analyses may use generalized additive or GAMLSS-based models to estimate nonlinear age effects, sex differences, peak ages, and cohort-adjusted trajectories across childhood, adulthood and aging.

Prediction/

Scripts for predictive modelling using CBSS features.

This folder may include phenotype prediction, cognitive prediction, or other machine-learning analyses based on fiber-wise, network-wise, or PCA-derived CBSS representations.

SC_FC/

Scripts for structure-function coupling analyses.

This folder is used to evaluate how CBSS-derived structural connectivity features relate to functional connectivity patterns, including SC-to-FC prediction or network-pair correspondence analyses.

Visualization/

Scripts for generating figures, atlas visualizations, network-level plots, fiber-profile plots, and other graphical summaries.

Typical Workflow

A typical CBSS workflow proceeds as follows:

1. Prepare atlas and template files.
2. Generate local neighborhood coordinates around atlas streamline points.
3. Project subject-specific diffusion maps onto atlas streamlines.
4. Aggregate projected values into fiber-wise or network-wise profiles.
5. Elastically align subject profiles to atlas mean profiles.
6. Derive final CBSS representations, such as mean, maximum, or PCA features.
7. Use the representations for lifespan, aging, prediction, SC-FC and visualization analyses.
Dependencies

The pipeline uses a combination of Python, R, MATLAB and C/MEX code.

Commonly used software and libraries include:

Python: numpy, nibabel, dipy, h5py, plotly
R: data.table, R.matlab, gamlss, gamlss.dist, gamlss.add, itsadug, pls, glmnet
MATLAB: batch-processing scripts, elastic-registration workflow orchestration, and C/MEX dynamic-programming functions
Tractography: TractFlow Singularity containers
Genetic analyses, where applicable: GCTA and fastGWA
Notes
The template/ files define the common CBSS atlas space and should remain unchanged unless a new atlas version is being generated.
The pipeline/ scripts should be run sequentially from step0 to step4.
The Analysis/ scripts assume that the final CBSS representations have already been generated.
File paths in the scripts may need to be updated before running on a new computing environment.

A shorter repository description could be:

```markdown
This repository provides the CBSS processing and analysis code. The `pipeline/` folder generates
