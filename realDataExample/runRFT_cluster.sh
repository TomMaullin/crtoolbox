#!/bin/bash

# -----------------------------------------------------------------------
# Work out path to repository
# -----------------------------------------------------------------------
RealPath() {
    (echo $(cd $(dirname "$1") && pwd -P)/$(basename "$1"))
}

CONFSETS_PATH=$(dirname $(RealPath "${BASH_SOURCE[0]}"))

# -----------------------------------------------------------------------
# Get number of slices
# -----------------------------------------------------------------------
nSlices=$(fslinfo /well/nichols/users/inf852/RFT_Ttest/COPE_diff_BODY_123117.nii | grep -m 1 dim3 | awk '{print $2}')
echo $nSlices

# Loop through all slices
sliceNo=0
while [ "$sliceNo" -lt "$nSlices" ]; do


	# -----------------------------------------------------------------------
	# Run slice jobs
	# -----------------------------------------------------------------------
	echo "Running slice "$sliceNo"..."

	# Submit config generation job
	fsl_sub -l $CONFSETS_PATH/realDataExample/results/log/ -N slice$sliceNo bash $CONFSETS_PATH/runRFT.sh $sliceNo > \
	/tmp/$$ && sliceID=$(awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}' /tmp/$$)
  
  	sliceNo=$(($sliceNo + 1))

done

echo "Job submitted, please monitor progress using qstat."
