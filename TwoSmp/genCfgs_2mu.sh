#!/bin/bash
fslpython -c "from TwoSmp.genCfgs_2mu import *; generateCfgs('$1',$2)"