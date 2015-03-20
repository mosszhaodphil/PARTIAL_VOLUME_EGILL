#!/bin/sh

# gen_pv_map anat T1t

file_name_perfusion_gm=$1
file_name_M0t_gm=$2

file_name_perfusion_scaled=perfusion_scaled.nii.gz

echo "Scaling..."

fslmaths $file_name_perfusion_gm -div $file_name_M0t_gm -mul 0.9 -mul 6000 $file_name_perfusion_scaled

echo Results saved in $file_name_perfusion_scaled

echo "Complete"