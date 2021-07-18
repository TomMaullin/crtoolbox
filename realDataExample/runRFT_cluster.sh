#!/bin/bash

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
# Run a slice
# -----------------------------------------------------------------------
echo "Running slice..."

# Submit config generation job
fsl_sub -l $CONFSETS_PATH/realDataExample/results/log/ -N slice bash $CONFSETS_PATH/runRFT.sh > \
/tmp/$$ && sliceID=$(awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}' /tmp/$$)

echo "Job submitted, please monitor progress using qstat."

