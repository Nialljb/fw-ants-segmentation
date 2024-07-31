#!/bin/sh

```
# Module: For running the ANTs pipeline for segmenting infant brain images on Flywheel
# Author: Chiara Casella, Niall Bourke


# Overview:
# This script is designed to run the ANTs pipeline for segmenting infant brain images on Flywheel. The pipeline consists of the following steps:
# 1. Register the input image to a template image
# 2. Register age-specific tissue and CSF priors from the BCP atlas to the input image via the template. The tissue prior was obtained by summing white matter and grey matter priors.
# 3. Segment the input image in native space using ANTs Atropos, two classes.
# 4. Use FSL to segment ventricles and subcortical grey matter. Rough masks were first drawn in template space and then moved to native space.
# 5. Extract volume estimates from the segmentations


# Usage:
# This script is designed to be run as a Flywheel Gear. The script takes two inputs:
# 1. The input image to segment
# 2. The age of the template to use in months (e.g. 3, 6, 12, 24, 48, 72)

# The script assumes that the input image is in NIfTI format. The script outputs the segmentations in native space.


# NOTE:
# Need txt output of volumes
# clean up intermediate files
# slicer bet & segmentations in native space (-A)

```

# Initialise the FSL environment
. ${FSLDIR}/etc/fslconf/fsl.sh

# Add c3d to the path
export PATH=$PATH:/flywheel/v0/utils/c3d-1.1.0/bin/c3d_affine_tool

#Define inputs
input_file=$1
age=$2

# Define the paths
FLYWHEEL_BASE=/flywheel/v0
INPUT_DIR=$FLYWHEEL_BASE/input/
WORK_DIR=$FLYWHEEL_BASE/work
OUTPUT_DIR=$FLYWHEEL_BASE/output
TEMPLATE_DIR=$FLYWHEEL_BASE/app/templates/${age}/
CONTAINER='[flywheel/ants-segmentation]'
template=${TEMPLATE_DIR}/template_${age}_degibbs.nii.gz

echo "permissions"
ls -ltra /flywheel/v0/

##############################################################################
# Handle INPUT file
# Check that input file exists

if [[ -e $input_file ]]; then
  echo "${CONTAINER}  Input file found: ${input_file}"

    # Determine the type of the input file
  if [[ "$input_file" == *.nii ]]; then
    type=".nii"
  elif [[ "$input_file" == *.nii.gz ]]; then
    type=".nii.gz"
  fi
  # Get the base filename
  # base_filename=`basename "$input_file" $type`
  native_img=`basename $input_file`
  native_img="${native_img%.nii.gz}"

else
  echo "${CONTAINER} no inputs were found within input directory $INPUT_DIR"
  exit 1
fi

##############################################################################

echo -e "\n --- Step 1: Register image to template --- "

# Define outputs in the following steps
native_bet_image=${OUTPUT_DIR}/native_bet_image.nii.gz
native_brain_mask=${OUTPUT_DIR}/native_brain_mask.nii.gz
template_brain_mask=${TEMPLATE_DIR}/brainMask.nii.gz

#bet image to help with registration to template
mri_synthstrip -i ${input_file} -o ${native_bet_image} -m ${native_brain_mask}
#bet ${input_file} ${native_bet_image}
#fslmaths ${native_bet_image} -bin ${native_brain_mask}
echo "BET image and mask created with synthstrip"
ls ${native_bet_image} ${native_brain_mask}
echo "***"  

# Register native BET image to template brain
echo -e "\n Run SyN registration"
#antsRegistrationSyN.sh -d 3 -i ${WORK_DIR}/itk.txt -t 'so' -f ${template} -m ${native_bet_image} -j 1 -o ${WORK_DIR}/bet_ -n 6
antsRegistrationSyN.sh -d 3 -t 's' -f ${template} -m ${native_bet_image} -j 1 -o ${WORK_DIR}/bet_ -n 6
echo "antsRegistrationSyN done"

# Short pause of 3 seconds
sleep 3

echo -e "\n --- Step 2: Apply registration to non-betted image --- "
# Get the affine and warp files from the registration
AFFINE_TRANSFORM=$(ls ${WORK_DIR}/*0GenericAffine.mat)
WARP=$(ls ${WORK_DIR}/*1Warp.nii.gz)
INVERSE_WARP=$(ls ${WORK_DIR}/*1InverseWarp.nii.gz)


# Transform priors (template space) to each subject's native space
echo "Transforming priors to native space for segmentation"

items=(
    "${TEMPLATE_DIR}/prior1.nii.gz"
    "${TEMPLATE_DIR}/prior2.nii.gz"
)

for item in "${items[@]}"; do
item_name=$(basename "$item" .nii.gz)
output_prior="${item_name}.nii.gz"
antsApplyTransforms -d 3 -i "${item}" -r ${native_bet_image} -o ${WORK_DIR}/"${output_prior}" -t ["$AFFINE_TRANSFORM",1] -t "${INVERSE_WARP}"
echo "$item_name transformed and saved to ${output_prior}"
done

# Short pause of 3 seconds
sleep 3


# Transform template mask and ventricles and subcortical grey matter masks (template space) to each subject's native space
echo "Transforming masks to native space"

items=(
    "${TEMPLATE_DIR}/ventricles_mask.nii.gz"
    "${TEMPLATE_DIR}/sub_GM_mask.nii.gz"
)

for item in "${items[@]}"; do
item_name=$(basename "$item" .nii.gz)
output_mask="${item_name}.nii.gz"
antsApplyTransforms -d 3 -i "${item}" -r ${native_bet_image} -o ${WORK_DIR}/"${output_mask}" -n NearestNeighbor -t ["$AFFINE_TRANSFORM",1] -t "${INVERSE_WARP}"
echo "$item_name transformed and saved to ${output_mask}"
done

# Short pause of 3 seconds
sleep 3


echo -e "\n--- Step 3: Segment image in native space with antsAtropos (2 priors) ---" 

#dilate native brain mask
fslmaths ${native_brain_mask} -dilM ${WORK_DIR}/native_brain_mask_dil.nii.gz

# Run Atropos
antsAtroposN4.sh -d 3 -a ${input_file} -x ${WORK_DIR}/native_brain_mask_dil.nii.gz -p ${WORK_DIR}/prior%d.nii.gz -c 2 -y 1 -y 2 -w 0.5 -o ${WORK_DIR}/ants_atropos_ -z 1

#get rid of extra tissue
fslmaths ${WORK_DIR}/ants_atropos_SegmentationPosteriors1.nii.gz -mul ${native_brain_mask} ${WORK_DIR}/ants_atropos_SegmentationPosteriors1_corrected.nii.gz   
fslmerge -t ${OUTPUT_DIR}/merged_priors.nii.gz ${WORK_DIR}/zero_filled_image.nii.gz ${WORK_DIR}/ants_atropos_SegmentationPosteriors1_corrected.nii.gz ${WORK_DIR}/ants_atropos_SegmentationPosteriors2.nii.gz
fslmaths ${OUTPUT_DIR}/merged_priors.nii.gz -Tmaxn ${OUTPUT_DIR}/${native_img}_atlas_2classes.nii.gz


echo -e "\n --- Step 4: Extract ventricles and subcortical grey matter --- "

# # Define posterior images from Atropos segmentation (segmentation in native space with 2 priors)
Posterior1=${WORK_DIR}/ants_atropos_SegmentationPosteriors1_corrected.nii.gz
Posterior2=${WORK_DIR}/ants_atropos_SegmentationPosteriors2.nii.gz

#ventricles
fslmaths ${WORK_DIR}/ventricles_mask.nii.gz -ero -dilM ${WORK_DIR}/refined_ventricles_mask.nii.gz
fslmaths ${Posterior2} -mul ${WORK_DIR}/refined_ventricles_mask.nii.gz ${WORK_DIR}/ventricles_mask_mul
fslmerge -t ${WORK_DIR}/merged_priors.nii.gz ${Posterior1} ${Posterior2} ${WORK_DIR}/ventricles_mask_mul.nii.gz
fslmaths ${WORK_DIR}/merged_priors.nii.gz -Tmean -mul $(fslval ${WORK_DIR}/merged_priors.nii.gz dim4) ${WORK_DIR}/merged_priors_Tsum
fslmaths ${WORK_DIR}/merged_priors_Tsum.nii.gz -thr 1.1 -bin ${WORK_DIR}/subtractmask



fslmaths ${Posterior2} -mul ${WORK_DIR}/subtractmask.nii.gz ${WORK_DIR}/ventricles_prior
fslmaths ${Posterior2} -sub ${WORK_DIR}/ventricles_prior.nii.gz ${WORK_DIR}/csf_prior

#create atlas
fslmaths ${native_brain_mask} -mul 0 ${WORK_DIR}/zero_filled_image.nii.gz
fslmerge -t ${WORK_DIR}/merged_2classes_priors.nii.gz ${WORK_DIR}/zero_filled_image.nii.gz ${WORK_DIR}/csf_prior.nii.gz ${WORK_DIR}/ventricles_prior.nii.gz ${Posterior2}
fslmaths ${WORK_DIR}/merged_2classes_priors.nii.gz -Tmaxn ${WORK_DIR}/atlas_2classes.nii.gz

#subcortical GM
fslmaths "${Posterior1}" -mul ${WORK_DIR}/sub_GM_mask.nii.gz ${WORK_DIR}/sub_GM_mask_mul
fslmerge -t ${WORK_DIR}/merged_sub_priors.nii.gz ${Posterior1} ${WORK_DIR}/csf_prior.nii.gz ${WORK_DIR}/ventricles_prior.nii.gz ${WORK_DIR}/sub_GM_mask_mul.nii.gz
fslmaths ${WORK_DIR}/merged_sub_priors.nii.gz -Tmean -mul $(fslval ${WORK_DIR}/merged_sub_priors.nii.gz dim4) ${WORK_DIR}/merged_sub_priors_Tsum
fslmaths ${WORK_DIR}/merged_sub_priors_Tsum.nii.gz -thr 1.1 -bin ${WORK_DIR}/subtractmask
fslmaths ${Posterior1} -mul ${WORK_DIR}/subtractmask.nii.gz ${WORK_DIR}/sub_GM_prior
fslmaths ${Posterior1} -sub ${WORK_DIR}/sub_GM_prior.nii.gz ${WORK_DIR}/Tissue_non_sub_GM_prior

#create atlas
fslmaths ${WORK_DIR}/brainMask_dil.nii.gz -mul 0 ${WORK_DIR}/zero_filled_image.nii.gz
fslmerge -t ${WORK_DIR}/merged_4classes_priors.nii.gz ${WORK_DIR}/zero_filled_image.nii.gz ${WORK_DIR}/Tissue_non_sub_GM_prior.nii.gz ${WORK_DIR}/sub_GM_prior.nii.gz ${WORK_DIR}/csf_prior.nii.gz ${WORK_DIR}/ventricles_prior.nii.gz
fslmaths ${WORK_DIR}/merged_4classes_priors.nii.gz -Tmaxn ${OUTPUT_DIR}/${native_img}_atlas_4classes.nii.gz


echo -e "\n --- Step 5: Run slicer and extract volume estimation from segmentations --- "
slicer ${native_bet_image} ${native_bet_image} -a ${WORK_DIR}/slicer_bet.png
slicer ${WORK_DIR}/ants_atropos_SegmentationPosteriors1_corrected.nii.gz ${WORK_DIR}/ants_atropos_SegmentationPosteriors1_corrected.nii.gz -a ${WORK_DIR}/slicer_seg1.png
slicer ${WORK_DIR}/ants_atropos_SegmentationPosteriors2.nii.gz ${WORK_DIR}/ants_atropos_SegmentationPosteriors2.nii.gz -a ${WORK_DIR}/slicer_seg2.png
pngappend ${WORK_DIR}/slicer_bet.png - ${WORK_DIR}/slicer_seg1.png - ${WORK_DIR}/slicer_seg2.png ${OUTPUT_DIR}/qc_montage.png

slicer ${WORK_DIR}/ants_atropos_SegmentationPosteriors1.nii.gz ${OUTPUT_DIR}/${native_img}_atlas_2classes.nii.gz -a ${OUTPUT_DIR}/slicer_atlas_2classes.png
slicer ${WORK_DIR}/ants_atropos_SegmentationPosteriors1.nii.gz ${OUTPUT_DIR}/${native_img}_atlas_4classes.nii.gz -a ${OUTPUT_DIR}/slicer_atlas_4classes.png


# Extract volumes of segmentations before extracting ventricles and subcortical GM
atlas=${OUTPUT_DIR}/${native_img}_atlas_2classes.nii.gz
declare -A labels=( ["tissue"]=1 ["csf"]=2)
output_csv="Tissue_and_csf_volumes.csv"
echo "Volume" > ${OUTPUT_DIR}/"$output_csv"
for tissue in "${!labels[@]}"; do
    label=${labels[$tissue]}
    mask_file=${OUTPUT_DIR}/"${tissue}_mask.nii.gz"
    fslmaths "$atlas" -thr $label -uthr $label -bin "$mask_file"
    volume=$(fslstats $atlas -k "$mask_file" -V | awk '{ print $2 }')  
    echo "${tissue},${volume}" >> ${OUTPUT_DIR}/"$output_csv"
done


# Extract volumes of segmentations after extracting ventricles and subcortical GM
atlas=${OUTPUT_DIR}/${native_img}_atlas_4classes.nii.gz
declare -A labels=( ["tissue_non_sub_GM"]=1 ["sub_GM"]=2 ["csf_non_ventricles"]=3 ["ventricles"]=4)
output_csv="All_volumes.csv"
echo "Volume" > ${OUTPUT_DIR}/"$output_csv"
for tissue in "${!labels[@]}"; do
    label=${labels[$tissue]}
    mask_file=${OUTPUT_DIR}/"${tissue}_mask.nii.gz"
    fslmaths "$atlas" -thr $label -uthr $label -bin "$mask_file"
    volume=$(fslstats $atlas -k "$mask_file" -V | awk '{ print $2 }')  
    echo "${tissue},${volume}" >> ${OUTPUT_DIR}/"$output_csv"
done


# --- Handle exit status --- #
# Check if the output directory is empty
if [ -z "$(find "$OUTPUT_DIR" -mindepth 1 -print -quit 2>/dev/null)" ]; then
    echo "Error: Output directory is empty"
    exit 1
fi







