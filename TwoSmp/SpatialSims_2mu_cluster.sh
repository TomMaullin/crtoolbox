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
fsl_sub -l $config_logdir -N setup bash $CONFSETS_PATH/TwoSmp/genCfgs_2mu.sh \
$CONFSETS_PATH/results $simNo > \
/tmp/$$ && cfgGenID=$(awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}' /tmp/$$)

# Wait until the baseline configuration file has been created
while [ ! -f $CONFSETS_PATH/results/sim$simNo/cfgs/baseline_cfg.yml ]; do sleep 1; done

echo "Configurations generated"

# -----------------------------------------------------------------------
# Run all configuration files
# -----------------------------------------------------------------------

echo "Submitting simulation instances..."

for cfg in $CONFSETS_PATH/results/sim$simNo/cfgs/cfg*.yml; do

	# Submit the job
	fsl_sub -j $cfgGenID -l $config_logdir -N setup \
	bash $CONFSETS_PATH/TwoSmp/SpatialSims_2mu.sh $CONFSETS_PATH/results/sim$simNo/cfgs/$cfg \
	> /tmp/$$ && setupID=$(awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}' /tmp/$$)

done

echo "Simulation instances submitted!"
