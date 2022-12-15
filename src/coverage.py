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


def check_violations(Fc, FcHat_plus, FcHat_minus, muHats, mus, sigmas, c, tau, a):

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
    binary_success = 1-(np.any(FcHatp_sub_Fc,axis=(1,2)) | np.any(Fc_sub_FcHatm,axis=(1,2)))# : AXES WONT WORK FOR 3D ATM

    # -------------------------------------------------------------------
    # Get stat along the Fc boundary
    # -------------------------------------------------------------------

    # Obtain statistic, g
    g = ((muHats-c)/(sigmas*tau))

    # -------------------------------------------------------------------
    # Boundary locations for FcHat
    # -------------------------------------------------------------------

    # Get boolean map for the boundary of Fc
    Fc_bdry_map = get_bdry_maps(mu, c)

    # Get coordinates for the boundary of Fc
    Fc_bdry_locs = get_bdry_locs(Fc_bdry_map)

    # Empty dict to store g along true boundary
    g_dFc = {}

    # Loop through fields
    for i in (np.arange(m)+1):

        # Get the values for gi along dFc
        g_dFc[str(i)] = get_bdry_values_concat(g[i,...], Fc_bdry_locs)
        
        # -------------------------------------------------------------------
        # Mu along dFc
        # -------------------------------------------------------------------

        # Obtain Mu along Fc
        mu_dFc[str(i)] = get_bdry_values_concat(mus[i,...], Fc_bdry_locs)

    # Empty dict to store interpolate g along true boundary
    g_dFc_interp = {}


    # -------------------------------------------------------------------
    # Boundary partitions 
    # -------------------------------------------------------------------

    # Get list of possible alphas to be considered
    alphas=list(powerset(np.arange(m)+1))

    # Locations of the dalpha boundaries
    dalphaFc_locs = {}
    dalphaFcHat_locs = {}

    # Loop through alphas
    for alpha in alphas:

        # Loop through possible elements of alpha
        for i in (np.arange(m)+1):

            # Get indices representing where mui and muhati are less
            # than or equal to c if i is in alpha.
            if i in alpha:

                # Mu
                in_dalphaFc_i = (mu_dFc[str(i)][:,1] <= c)

            # Get indices representing where mui is greater than c
            # if i is not in alpha.
            else:

                # Mu
                in_dalphaFc_i = (mu_dFc[str(i)][:,1] > c)

            # Get a boolean index telling us whether each location in dFc
            # belongs to d^\alpha Fc.
            if i == 1:

                # Initial running product for mu
                in_dalphaFc = np.array(in_dalphaFc_i)

            else:

                # Update running product for mu
                in_dalphaFc = in_dalphaFc*in_dalphaFc_i

        # Convert boolean 1/0s into coordinates
        dalphaFc_loc = np.where(in_dalphaFc)[0]

        # Save locations
        dalphaFc_locs[np.array2string(alpha)] = dalphaFc_loc


    # -------------------------------------------------------------------
    # Get weights from mu and muhat along boundary partitions for 
    # interpolation
    # -------------------------------------------------------------------

    # Empty dicts for weights
    weights_dFc = {}

    # Loop through boundary partitions
    for alpha in alphas:

        # New empty dicts for weights
        weights_dalphaFc = {}

        # Loop through i in alpha getting values for interpolation
        for i in alpha:

            # Get mui and muhati along dalpha Fc
            mui_dalphaFc = mu_dFc_partitioned[np.array2string(alpha)][str(i)]

            # For true boundary
            dalphaFc_mui_weights = get_bdry_weights_concat(mui_dalphaFc, c)

            # Save weights
            weights_dalphaFc[str(i)] = dalphaFc_mui_weights

        # Save weights
        weights_dFc[np.array2string(alpha)] = weights_dalphaFc


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

    # -------------------------------------------------------------------
    # Work out whether we observed successful sets.
    # -------------------------------------------------------------------
    # Record if we saw a violation via the interpolation method
    interp_success = np.all(bdry_lowerCheck,axis=(1)) & np.all(bdry_upperCheck,axis=(1))

    # -------------------------------------------------------------------
    # Success as a whole
    # -------------------------------------------------------------------
    # Final success
    success = interp_success & binary_success

    # Return result
    return(success, {'Interp: ': interp_success, 'Binary: ': binary_success})
