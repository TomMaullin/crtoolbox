#!/bin/bash
#fslpython -c "from SpatialSims import *; SpatialSims('/users/nichols/inf852/tmp', $1, {'type': 'circle2D', 'center': np.array([0,0]), 'fwhm': np.array([3,3]), 'r': 30, 'mag': 3}, $2, 2, np.linspace(0,1,21))"
fslpython -c "from SpatialSims import *; SpatialSims('/users/nichols/inf852/tmp', $1, {'type': 'ramp2D', 'a': 1, 'b': 3, 'orient': 'horizontal'}, $2, 2, np.linspace(0,1,21))"