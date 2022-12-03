import os
import sys
import time
import numpy as np
from boundary import *
from fileio import *
from bootstrap import *
from setTheory2mu import *
import yaml
import matplotlib.pyplot as plt


def check_violations(Fc, FcHat_plus, FcHat_minus, muHats, sigmas, c, tau, a):

    # -------------------------------------------------------------------
    # Some set logic to work out violations
    # -------------------------------------------------------------------

    # Obtain FcHat^+\FcHat based on the estimated boundary. This variable
    # has axes corresponding to [pvalue, field dimensions]
    FcHatp_sub_Fc = FcHat_plus & ~Fc

    # Obtain Fc\FcHat^- based on the estimated boundary. This variable
    # has axes corresponding to [pvalue, field dimensions]
    Fc_sub_FcHatm = Fc & ~FcHat_minus

    # Record if we saw a violation in the estimated boundary based sets
    estBdry_success = 1-(np.any(FcHatp_sub_Fc,axis=(1,2)) | np.any(Fc_sub_FcHatm,axis=(1,2)))# : AXES WONT WORK FOR 3D ATM

    # -------------------------------------------------------------------
    # Get stat along the Fc boundary
    # -------------------------------------------------------------------

    # Empty dict to store g along true boundary
    g_dFc = {}

    # Loop through fields
    for i in np.arange(m):

        # Get the values for gi along dFc
        g_dFc[str(i+1)] = get_bdry_values_concat(g[i,...], Fc_bdry_locs)

    # Empty dict to store interpolate g along true boundary
    g_dFc_interp = {}

    # Boolean to tell us if this is the first alpha we've looked at
    Firstalpha = True

    # Interpolate g
    for alpha in alphas:

        # Initial structure to hold interpolated g
        g_dalphaFc_interp = {}

        # Boolean to tell us if this is the first value of i we're looking at
        Firsti = True

        # Check if we have recorded values for this boundary
        if np.array2string(alpha) in dalphaFc_locs:

            # Loop through all fields relevent to boundary.
            for i in alpha:

                # Get gi on dFc
                gi_dFc = g_dFc[str(i)] 

                # Get locations for dalphaFc
                dalphaFc_loc = dalphaFc_locs[np.array2string(alpha)]

                # Get gi on dalphaFc
                gi_dalphaFc = gi_dFc[dalphaFc_loc,:]

                # Get weights for interpolation
                dalphaFc_weightsi = weights_dFc[np.array2string(alpha)][str(i)]

                # Interpolate gi
                gi_dalphaFc_interp = get_bdry_vals_interpolated_concat(gi_dalphaFc,dalphaFc_weightsi)

                # If this is the first i we've looked at in alpha
                if Firsti:

                    # set the minimum interpolated values to gi_dalphaFc_interp
                    ming_dalphaFc_interp = gi_dalphaFc_interp 

                    # Set first to false
                    Firsti = False

                else:

                    # Take the minimum between current gi_dalphaFc_interp and previous
                    ming_dalphaFc_interp = np.minimum(ming_dalphaFc_interp, gi_dalphaFc_interp) 

            # -------------------------------------------------------------------
            # Get stat along the Fc boundary
            # -------------------------------------------------------------------

            # If this is the first alpha we've looked at, set the stat array to ming_dalphaFc_interp
            if Firstalpha:

                stat_dFc = ming_dalphaFc_interp

                # Now we've seen an alpha
                Firstalpha = False

            # Otherwise concatenate it on
            else:

                stat_dFc = np.concatenate((ming_dalphaFc_interp,stat_dFc),axis=-1)

    # -------------------------------------------------------------------
    # Check whether there were any boundary violations using interpolated
    # boundary values (checking if voxels had values corresponding to no
    # violations, etc)
    # -------------------------------------------------------------------

    # Perform lower check on stat map using thresholds based on the
    # estimated boundary
    bdry_lowerCheck = stat_dFc >= a[:,0,:,0]

    # Perform upper check on stat map using thresholds based on the
    # estimated boundary
    bdry_upperCheck = stat_dFc <= a[:,1,:,0]

    # Perform lower check on stat map using thresholds based on the
    # true boundary
    bdry_lowerCheck_trueBdry = stat_dFc >= a_trueBdry[:,0,:,0]

    # Perform upper check on stat map using thresholds based on the
    # true boundary
    bdry_upperCheck_trueBdry = stat_dFc <= a_trueBdry[:,1,:,0]
