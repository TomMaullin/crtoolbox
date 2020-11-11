#!/bin/bash
fslpython -c "from SpatialSims import *; SpatialSims('/users/nichols/inf852/tmp', $1, 'rampHoriz2D', $2, 2, np.linspace(0,1,20))"