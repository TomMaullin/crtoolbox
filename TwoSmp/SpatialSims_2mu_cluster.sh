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

# include parse_yaml function
. $CONFSETS_PATH/parse_yaml.sh

# -----------------------------------------------------------------------
# Generate configuration files
# -----------------------------------------------------------------------
echo "Generating configurations..."

# Submit config generation job
fsl_sub -l $CONFSETS_PATH/results/sim$simNo/log/ -N genCfgs bash $CONFSETS_PATH/genCfgs_2mu.sh \
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

	# read yaml file to get output directory
	eval $(parse_yaml $cfg "config_")

	# Submit the job
	fsl_sub -j $cfgGenID -l $CONFSETS_PATH/results/sim$simNo/log/ -N cfg$config_cfgId \
	bash $CONFSETS_PATH/SpatialSims_2mu.sh $cfg \
	> /tmp/$$ && simInstancesID=$(awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}' /tmp/$$),$simInstancesID

done

# -----------------------------------------------------------------------
# Run concatenation job
# -----------------------------------------------------------------------

echo "Submitting concatenation job..."

# Submit the job
fsl_sub -j $simInstancesID -l $CONFSETS_PATH/results/sim$simNo/log/ -N concat \
bash $CONFSETS_PATH/join_and_plot_2mu.sh $CONFSETS_PATH/results $simNo \
> /tmp/$$ && concatID=$(awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}' /tmp/$$)

echo "All jobs submitted, please monitor progress using qstat."

