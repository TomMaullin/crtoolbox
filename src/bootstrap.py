import os
import sys
import time
import numpy as np
from boundary import *
from setTheory2mu import *
from fileio import *
import yaml
import matplotlib.pyplot as plt
# m, nReals, pVals, resids (resids_dFcHat_partitioned)

def bootstrap_resids(resids, m, n_boot, p, n_sub, interp_weights):

    # Work out number of p values
    nPvals = len(p)

    # Get list of possible alphas to be considered
    alphas=list(powerset(np.arange(m)+1))

    # Dimensions of bootstrap variables
    boot_dim = np.array([n_sub, 1]) 

    # -------------------------------------------------------------------
    # Bootstrap 
    # -------------------------------------------------------------------
    # Initialize empty bootstrap stores for supremum along each boundary
    # segment
    min_supg_dFcHat = {}

    # Loop through boundary segments
    for alpha in alphas:

        # Save empty bootstrap stores
        min_supg_dFcHat[np.array2string(alpha)] = np.zeros(n_boot)

    # Save empty bootstrap stores
    min_supg_dFcHat['max'] = np.zeros(n_boot)

    # For each bootstrap record the max of the residuals along the
    # boundary
    for b in np.arange(n_boot):

        # Obtain bootstrap variables
        boot_vars = 2*np.random.randint(0,2,boot_dim)-1

        # Reshape for broadcasting purposes (extra axis refers to the fact we have
        # inner and outer boundary values in the last axes of resid_FsdGc_bdry_concat
        # and resid_FsdGcHat_bdry_concat)
        boot_vars = boot_vars.reshape((*boot_vars.shape),1)
        
        # -------------------------------------------------------------------------
        # Bootstrap residuals along dalphaFc
        # -------------------------------------------------------------------------
        boot_resids_dFcHat = {}

        # Loop through boundary partitions
        for alpha in alphas:

            # New empty dicts for FcHat
            boot_resids_dalphaFcHat = {}

            # New empty dict for gi along dalphaFcHat
            boot_g_dalphaFcHat = {}

            # Loop through i in alpha getting gi interpolated
            for i in alpha:

                # ------------------------------------------------------
                # Get residuals true and estimated boundary
                # ------------------------------------------------------

                # Get original residuals of mui along boundary Aci
                residsi_dalphaFcHat = resids[np.array2string(alpha)][str(i)]

                # ------------------------------------------------------
                # Bootstrap residuals
                # ------------------------------------------------------

                # Multiply by rademacher variables
                boot_residsi_dalphaFcHat = boot_vars*residsi_dalphaFcHat

                # ------------------------------------------------------
                # Get gi along dalpha FcHat
                # ------------------------------------------------------

                # Sum across subjects to get the bootstrapped a values along
                # the boundary of dalphaFcHat. (Note: For some reason this is 
                # much faster if performed seperately for each of the last rows. 
                # I am still looking into why this is)
                boot_gi_dalphaFcHat = np.zeros(boot_residsi_dalphaFcHat.shape[-2:])
                boot_gi_dalphaFcHat[...,0] = np.sum(boot_residsi_dalphaFcHat[...,0], axis=0)/np.sqrt(n_sub)
                boot_gi_dalphaFcHat[...,1] = np.sum(boot_residsi_dalphaFcHat[...,1], axis=0)/np.sqrt(n_sub)

                # Obtain bootstrap standard deviations for muHat i along dalphaFcHat. 
                # (Note: For some reason this is much faster if performed seperately
                # for each of the last rows. I am still looking into why this is)
                boot_sigmai_dalphaFcHat = np.zeros(boot_residsi_dalphaFcHat.shape[-2:])
                boot_sigmai_dalphaFcHat[...,0] = np.std(boot_residsi_dalphaFcHat[...,0], axis=0, ddof=1)
                boot_sigmai_dalphaFcHat[...,1] = np.std(boot_residsi_dalphaFcHat[...,1], axis=0, ddof=1)

                # Divide by the boostrap standard deviation of mui
                boot_gi_dalphaFcHat = boot_gi_dalphaFcHat/boot_sigmai_dalphaFcHat

                # ------------------------------------------------------
                # Interpolate along dalpha FcHat
                # ------------------------------------------------------

                # Get weights
                dalphaFcHat_muHati_bdry_weights = interp_weights[np.array2string(alpha)][str(i)]

                # Interpolation for gi along dalphaFc
                boot_gi_dalphaFcHat = get_bdry_vals_interpolated_concat(boot_gi_dalphaFcHat,dalphaFcHat_muHati_bdry_weights)

                # ------------------------------------------------------
                # Get supremum along dalphaFcHat
                # ------------------------------------------------------

                # Get the maximum of gi along dalphaFcHat if it exists
                if np.prod(boot_gi_dalphaFcHat.shape) > 0:

                    boot_g_dalphaFcHat[str(i)] = boot_gi_dalphaFcHat

            # -----------------------------------------------------------------
            # Get sup_{dalpha FcHat} |min_{alpha} g^i|
            # -----------------------------------------------------------------

            # Reset first counter
            first = True

            # Loop through i in alpha getting min(gi)
            for i in alpha:

                # Check if we have boundary values for gi
                if str(i) in boot_g_dalphaFcHat:

                    # If this is the first time initalize the min(g) array
                    if first:

                        # Initialize elementwise minimum across gi of boot_g
                        boot_ming_dalphaFcHat = boot_g_dalphaFcHat[str(i)] 

                        # We are no longer looking at the first
                        first = False

                    else:

                        # Get elementwise minimum across gi of boot_g
                        boot_ming_dalphaFcHat = np.minimum(boot_ming_dalphaFcHat,boot_g_dalphaFcHat[str(i)])

            # If we saw some boundary values (i.e. we saw at least the first i),
            # then the array was initialized
            if not first:

                # Save sup(|min(gi)|) along estimated dalpha FcHat
                boot_sup_ming_dalphaFcHat = np.max(np.abs(boot_ming_dalphaFcHat))

                # Save bootstrap result
                min_supg_dFcHat[np.array2string(alpha)][b] = boot_sup_ming_dalphaFcHat

        # -----------------------------------------------------------------
        # Get max_{P(M)} sup_{dalpha Fc} |min_{alpha} g^i| and
        # max_{P(M)} sup_{dalpha FcHat} |min_{alpha} g^i| 
        # -----------------------------------------------------------------

        # Loop through alphas taking maximum
        for alpha in alphas:

            # Check if we have boundary values for dalpha FcHat
            if np.array2string(alpha) in min_supg_dFcHat:

                # Update the minimum value we've seen
                min_supg_dFcHat['max'][b] = np.maximum(min_supg_dFcHat['max'][b],min_supg_dFcHat[np.array2string(alpha)][b])

    # -------------------------------------------------------------------
    # Obtaining a from percentiles of the max distribution
    # -------------------------------------------------------------------

    # Drop the instances where the boundary length was zero
    min_supg_dFcHat['max'] = min_supg_dFcHat['max'][min_supg_dFcHat['max']!=0]

    # If we have recorded values get their quantiles
    if (np.prod(min_supg_dFcHat['max'].shape) > 0):
        # Get the a estimates for the estimated boundary
        a = np.percentile(min_supg_dFcHat['max'], 100*p).reshape(nPvals,1,1,1) # [pvals, 1, [1 for _ in dim>1]]
    else:
        # Set to inf by default
        a = np.Inf*np.ones((nPvals,1,1,1))

    # Reformat them to an array form useful for boolean operation
    a = np.concatenate((-a,a),axis=1)

    return(a)

