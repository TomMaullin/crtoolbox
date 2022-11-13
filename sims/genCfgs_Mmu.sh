#!/bin/bash
fslpython -c "from MSmp.genCfgs_Mmu import *; generateCfgs('$1',$2)"