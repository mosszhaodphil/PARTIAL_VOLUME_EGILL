#!/bin/sh

# gen_pv_map anat T1t

file_name_structure=$1
file_name_t1=$2

file_name_extracted_brain=structure_brain
file_name_registration=str_to_t1
file_name_seg_csf=str_to_t1_pve_0
file_name_seg_gm=str_to_t1_pve_1
file_name_seg_wm=str_to_t1_pve_2
file_name_csf_raw=pvcsf_raw
file_name_gm_raw=pvgm_raw
file_name_wm_raw=pvwm_raw
file_name_csf=pvcsf
file_name_gm=pvgm
file_name_wm=pvwm
file_name_csf_mask=csf_mask
file_name_gm_mask=gm_mask
file_name_wm_mask=wm_mask

echo "Performing partial volume tissue segmentation..."

# Brain extraction
bet $file_name_structure $file_name_extracted_brain

# Rigid registration
flirt -dof 6 -in $file_name_extracted_brain -ref $file_name_t1 -out $file_name_registration

# Tissue segmentation
fast -p $file_name_registration

# Resample
applywarp -i $file_name_seg_csf -r $file_name_t1 -o $file_name_csf_raw -s --interp=trilinear
applywarp -i $file_name_seg_gm -r $file_name_t1 -o $file_name_gm_raw -s --interp=trilinear
applywarp -i $file_name_seg_wm -r $file_name_t1 -o $file_name_wm_raw -s --interp=trilinear

# Thresholding at 0.1
fslmaths $file_name_csf_raw -thr 0.1 -min 1 $file_name_csf
fslmaths $file_name_gm_raw -thr 0.1 -min 1 $file_name_gm
fslmaths $file_name_wm_raw -thr 0.1 -min 1 $file_name_wm

# Binary tissue map
fslmaths $file_name_csf -bin $file_name_csf_mask
fslmaths $file_name_gm -bin $file_name_gm_mask
fslmaths $file_name_wm -bin $file_name_wm_mask


echo "Complete"