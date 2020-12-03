#!/bin/bash
#fslpython -c "from SpatialSims import *; SpatialSims('/users/nichols/inf852/tmp', $1, {'type': 'circle2D', 'center': np.array([0,0]), 'fwhm': np.array([3,3]), 'r': 30, 'mag': 3}, $2, 2, np.linspace(0,1,21))"
fslpython -c "from SpatialSims_2mu_v1 import *; SpatialSims_2mu_v1('/users/nichols/inf852/tmp', $1, $2, 2, np.linspace(0,1,21))"