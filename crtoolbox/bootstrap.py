import os
import sys
import time
import numpy as np
from crtoolbox.lib.boundary import *
from crtoolbox.lib.set_theory import powerset
from crtoolbox.lib.fileio import *
import yaml
import matplotlib.pyplot as plt

def bootstrap_resids(resid_vals, resid_weights, m, n_boot, p, n_sub):
    """
    Calculate bootstrap residuals along boundary segments and compute
    the quantile values based on the bootstrap distribution.
    
    Parameters:
    -----------
    resid_vals : dict
        Dictionary containing the residuals of boundary values.
    resid_weights : dict
        Dictionary containing the weights of residuals.
    m : int
        Number of conditions.
    n_boot : int
        Number of bootstrap iterations.
    p : array_like
        Percentile values for calculating quantiles.
    n_sub : int
        Number of subjects.
        
    Returns:
    --------
    a : ndarray
        Array of quantile values.
    """

    # Work out number of p values
    nPvals = len(p)

    # Get list of possible alphas to be considered
    alphas=list(powerset(np.arange(m)+1))

    # Dimensions of bootstrap variables
    boot_dim = np.array([n_sub, 1, 1]) 

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
        boot_vars = 2*np.random.randint(0,2,boot_dim,dtype="int8")-1
        
        # Loop through boundary partitions
        for alpha in alphas:
            
            # New empty dict for gi along dalphaFcHat
            boot_g_dalphaFcHat = {}

            # Loop through i in alpha getting gi interpolated
            for i in alpha:

                # ------------------------------------------------------
                # Get residuals true and estimated boundary
                # ------------------------------------------------------

                # Get original residuals of mui along boundary Aci
                residsi_dalphaFcHat = np.array(resid_vals[np.array2string(alpha)][str(i)], dtype="float32")

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
                mean_outer = np.sum(boot_residsi_dalphaFcHat[...,0], axis=0)/n_sub
                mean_inner = np.sum(boot_residsi_dalphaFcHat[...,1], axis=0)/n_sub

                # Get sum of squares
                ssq_outer = np.sum((boot_residsi_dalphaFcHat[...,0]-mean_outer)**2, axis=0)/(n_sub-1)
                ssq_inner = np.sum((boot_residsi_dalphaFcHat[...,1]-mean_inner)**2, axis=0)/(n_sub-1)

                # Obtain bootstrap g
                boot_gi_outer = np.sqrt(n_sub)*mean_outer/np.sqrt(ssq_outer)
                boot_gi_inner = np.sqrt(n_sub)*mean_inner/np.sqrt(ssq_inner)
                
                # clean up
                del mean_outer, mean_inner, ssq_outer, ssq_inner

                # ------------------------------------------------------
                # Interpolate along dalpha FcHat
                # ------------------------------------------------------

                # Get weights
                dalphaFcHat_muHati_bdry_weights = np.array(resid_weights[np.array2string(alpha)][str(i)], dtype="float32")

                # Work out weights
                outer_weights = dalphaFcHat_muHati_bdry_weights[...,0]
                inner_weights = dalphaFcHat_muHati_bdry_weights[...,1]
                

                # Work out interpolated values
                boot_gi_dalphaFcHat = inner_weights*boot_gi_inner + outer_weights*boot_gi_outer

                # ------------------------------------------------------
                # Get supremum along dalphaFcHat
                # ------------------------------------------------------

                # Get the maximum of gi along dalphaFcHat if it exists
                boot_g_dalphaFcHat[str(i)] = boot_gi_dalphaFcHat

                # clean up
                del boot_gi_dalphaFcHat

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
        a = np.percentile(min_supg_dFcHat['max'], 100*p).reshape(nPvals,1)
    else:
        # Set to inf by default
        a = np.Inf*np.ones((nPvals,1))

    # Reformat them to an array form useful for boolean operation
    a = np.concatenate((-a,a),axis=1)

    return(a)
