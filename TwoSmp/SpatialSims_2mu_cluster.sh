#!/bin/bash

# -----------------------------------------------------------------------
# Set simulation number
# -----------------------------------------------------------------------
simNo=1

# -----------------------------------------------------------------------
# Work out path to repository
# -----------------------------------------------------------------------
RealPath() {
    (echo $(cd $(dirname "$1") && pwd -P)/$(basename "$1"))
}

CONFSETS_PATH=$(dirname $(RealPath "${BASH_SOURCE[0]}"))

# -----------------------------------------------------------------------
# Generate configuration files
# -----------------------------------------------------------------------
echo "Generating configurations..."

# Submit config generation job
fsl_sub -l $CONFSETS_PATH/results/log/ -N genCfgs bash $CONFSETS_PATH/genCfgs_2mu.sh \
$CONFSETS_PATH/results $simNo > \
/tmp/$$ && cfgGenID=$(awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}' /tmp/$$)

# Wait until the baseline configuration file has been created
while [ ! -f $CONFSETS_PATH/results/sim$simNo/cfgs/baseline_cfg.yml ]; do sleep 1; done

echo "Configurations generated"

# -----------------------------------------------------------------------
# Run all configuration files
# -----------------------------------------------------------------------

echo "Submitting simulation instances..."

i=0
for cfg in $CONFSETS_PATH/results/sim$simNo/cfgs/cfg*.yml; do

	# Submit the job
	fsl_sub -j $cfgGenID -l $CONFSETS_PATH/results/log/ -N cfg$i \
	bash $CONFSETS_PATH/SpatialSims_2mu.sh $cfg \
	> /tmp/$$ && simInstancesID=$(awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}' /tmp/$$)

	i=$(($i + 1))

done

echo "Simulation instances submitted!"
