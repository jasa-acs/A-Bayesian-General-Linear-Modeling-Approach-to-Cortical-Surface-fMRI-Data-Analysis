cd ~/BayesianGLM/locations/

date

##################################################################
# SET UP RESAMPLING FROM 32K TO 6K (FOR BAYESIAN GLM)
##################################################################

# Step 1. Create new spherical GIFTI files using wb_command -surface-create-sphere  
~/workbench/bin_rh_linux64/wb_command -surface-create-sphere 6000 Sphere.6k.R.surf.gii
~/workbench/bin_rh_linux64/wb_command -surface-flip-lr Sphere.6k.R.surf.gii Sphere.6k.L.surf.gii
~/workbench/bin_rh_linux64/wb_command -set-structure Sphere.6k.R.surf.gii CORTEX_RIGHT
~/workbench/bin_rh_linux64/wb_command -set-structure Sphere.6k.L.surf.gii CORTEX_LEFT


# Step 2. Create CIFTI space with 6K vertices 

# a) Use -cifti-separate to get separate CORTEX_LEFT, CORTEX_RIGHT and both medial walls

# cifti.32K.dscalar.nii can be any dscalar file in the original 32K resolution (~32K locations per hemisphere)
date
~/workbench/bin_rh_linux64/wb_command -cifti-separate cifti.32K.dscalar.nii COLUMN -volume-all vol.32K.nii -label label.32K.nii -metric CORTEX_LEFT cifti.32K.L.shape.gii -roi roi.32K.L.shape.gii -metric CORTEX_RIGHT cifti.32K.R.shape.gii -roi roi.32K.R.shape.gii
date

# b) Use wb_commmand -metric-resample to create 6k versions of these using the 32k sphere and your new 6k sphere

~/workbench/bin_rh_linux64/wb_command -metric-resample cifti.32K.L.shape.gii Sphere.32K.L.surf.gii Sphere.6k.L.surf.gii BARYCENTRIC cifti.6K.L.shape.gii -current-roi roi.32K.L.shape.gii -valid-roi-out valid-roi.6K.L.shape.gii
~/workbench/bin_rh_linux64/wb_command -metric-resample cifti.32K.R.shape.gii Sphere.32K.R.surf.gii Sphere.6k.R.surf.gii BARYCENTRIC cifti.6K.R.shape.gii -current-roi roi.32K.R.shape.gii -valid-roi-out valid-roi.6K.R.shape.gii
~/workbench/bin_rh_linux64/wb_command -metric-resample roi.32K.L.shape.gii Sphere.32K.L.surf.gii Sphere.6k.L.surf.gii BARYCENTRIC roi.6K.L.shape.gii 
~/workbench/bin_rh_linux64/wb_command -metric-resample roi.32K.R.shape.gii Sphere.32K.R.surf.gii Sphere.6k.R.surf.gii BARYCENTRIC roi.6K.R.shape.gii 


# Step 3. Create a new CIFTI dense timeseries.

~/workbench/bin_rh_linux64/wb_command -cifti-create-dense-timeseries ts.6K.dtseries.nii -volume vol.32K.nii label.32K.nii -left-metric cifti.6K.L.shape.gii -roi-left roi.6K.L.shape.gii -right-metric cifti.6K.R.shape.gii -roi-right roi.6K.R.shape.gii
~/workbench/bin_rh_linux64/wb_command -cifti-create-dense-timeseries ts.32K.dtseries.nii -volume vol.32K.nii label.32K.nii -left-metric cifti.32K.L.shape.gii -roi-left roi.32K.L.shape.gii -right-metric cifti.32K.R.shape.gii -roi-right roi.32K.R.shape.gii


# Note: actual resampling performed in BayGLM1.m

date