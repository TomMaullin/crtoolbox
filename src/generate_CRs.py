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


# ===============================================================================
#
# This file contains the central functions for generating confidence regions for
# conjunction inference. This code was last updated on 13/11/2022
#
# Author: Tom Maullin 
# Contact: TomMaullin@gmail.com
#
# ===============================================================================


# ===============================================================================
#
#
#
# -------------------------------------------------------------------------------
#
# The below function takes the following inputs:
#
# -------------------------------------------------------------------------------
#
# data - shape (number of samples, number of  subjects, *dim)
#
# -------------------------------------------------------------------------------
#
# And returns the following outputs:
#
# ===============================================================================
def generate_CRs(data, c, p, mask=None, n_boot=5000, tau='1/np.sqrt(n_sub)'):


    # Work out m, the number of samples we are considering
    m = data.shape[0]

    # Get number of subjects
    n_sub = data.shape[1]

    # Get image dimensions
    image_dim = data.shape[2:]

    # Evaulate tau
    tau = eval(tau)

    # -------------------------------------------------------------------
    # Mean and standard deviation estimates
    # -------------------------------------------------------------------

    # Obtain mu estimate
    muHats = np.mean(data, axis=1).reshape(m,*image_dim)

    # Obtain sigma
    sigmas = np.std(data, axis=1).reshape(m,*image_dim)

    # -------------------------------------------------------------------
    # Boundary locations and values
    # -------------------------------------------------------------------

    # Make a structure to hold the estimated boundary weights in array form 
    est_bdry_weights_concat = {}

    # Make a structure to hold the true and estimated boundary locations
    est_bdry_locs = {}

    for i in np.arange(m):

        # Get muhat for this sample
        muHat = muHats[i,...]

        # -------------------------------------------------------------------
        # Boundary locations for AcHati
        # -------------------------------------------------------------------
        # Get boolean maps for the boundary of AcHat
        AcHat_bdry_map = get_bdry_maps(muHat, c)

        # Get coordinates for the boundary of AcHat
        AcHat_bdry_locs = get_bdry_locs(AcHat_bdry_map)

        # Save boundary locations
        est_bdry_locs['AcHat'+str(i+1)] = AcHat_bdry_locs

        # Delete map as we no longer need it
        del AcHat_bdry_map

        # -------------------------------------------------------------------
        # Interpolation weights for AcHati boundary (Array version)
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for AcHati
        AcHat_bdry_vals_concat = get_bdry_values_concat(muHat, AcHat_bdry_locs)

        # Obtain the weights along the boundary for AcHati
        AcHat_bdry_weights_concat = get_bdry_weights_concat(AcHat_bdry_vals_concat, c)

        # Save boundary weights
        est_bdry_weights_concat['AcHat'+str(i+1)] = AcHat_bdry_weights_concat

        # Delete values as we no longer need them
        del AcHat_bdry_vals_concat


    # -------------------------------------------------------------------
    # Get minimum fields
    # -------------------------------------------------------------------
    # This is named cap as the excursion set of the minimum field is
    # the intersection of all fields (\cap in latex)
    cap_muHat = np.amin(muHats,axis=0)

    # -------------------------------------------------------------------
    # Boundary locations for FcHat
    # -------------------------------------------------------------------
    # Get boolean map for the boundary of FcHat
    FcHat_bdry_map = get_bdry_maps(cap_muHat, c)

    # Get coordinates for the boundary of FcHat
    FcHat_bdry_locs = get_bdry_locs(FcHat_bdry_map)

    # Save boundary locations
    est_bdry_locs['FcHat'] = FcHat_bdry_locs

    # Delete maps as we no longer need them
    del FcHat_bdry_map

    # Empty dict to store residuals
    resids_dFcHat = {}

    # Empty dicts to store  muHat
    muHat_dFcHat = {}

    # Loop through to get residuals
    for i in np.arange(m):

        # -------------------------------------------------------------------
        # Residuals along dFcHat
        # -------------------------------------------------------------------

        # Obtain residuals
        resid = (data[i,...]-muHats[i,...])/sigmas[i,...]

        # Residuals along FcHat boundary
        resid_dFcHat_concat = get_bdry_values_concat(resid, FcHat_bdry_locs)

        # Save residuals
        resids_dFcHat['field'+str(i+1)] = resid_dFcHat_concat

        # -------------------------------------------------------------------
        # MuHat along dFcHat
        # -------------------------------------------------------------------

        # Obtain MuHat along FcHat
        muHat_dFcHat_concat = get_bdry_values_concat(muHats[i,...], FcHat_bdry_locs)

        # Save mu
        muHat_dFcHat['field'+str(i+1)] = muHat_dFcHat_concat

    # Delete data as it is longer needed
    del data, resid

    # -------------------------------------------------------------------
    # Boundary partitions 
    # -------------------------------------------------------------------

    # Get list of possible alphas to be considered
    alphas=list(powerset(np.arange(m)+1))

    # Locations of the dalpha boundaries
    dalphaFcHat_locs = {}

    # Loop through alphas
    for alpha in alphas:

        # Loop through possible elements of alpha
        for i in (np.arange(m)+1):

            # Get indices representing where muhati is less
            # than or equal to c if i is in alpha.
            if i in alpha:

                # Muhat
                in_dalphaFcHat_i = (muHat_dFcHat['field'+str(i)][:,1] <= c)

            # Get indices representing where muhati is greater than c
            # if i is not in alpha.
            else:

                # Muhat
                in_dalphaFcHat_i = (muHat_dFcHat['field'+str(i)][:,1] > c)


            # Get a boolean index telling us whether each location in dFc
            # belongs to d^\alpha Fc.
            if i == 1:

                # Initial running product for muHat
                in_dalphaFcHat = np.array(in_dalphaFcHat_i)

            else:

                # Update running product for muhat
                in_dalphaFcHat = in_dalphaFcHat*in_dalphaFcHat_i

        # Convert boolean 1/0s into coordinates
        dalphaFcHat_loc = np.where(in_dalphaFcHat)[0]

        # Save locations
        dalphaFcHat_locs[np.array2string(alpha)] = dalphaFcHat_loc

    # -------------------------------------------------------------------
    # Get residuals and muhat along boundary partitions
    # -------------------------------------------------------------------

    # Empty dicts for residuals and muHat
    resids_dFcHat_partitioned = {}
    muHat_dFcHat_partitioned = {}

    # Loop through boundary partitions
    for alpha in alphas:

        # New empty dicts for FcHat
        resids_dalphaFcHat = {}
        muHat_dalphaFcHat = {}

        # Get dalpha locations
        dalphaFcHat_loc = dalphaFcHat_locs[np.array2string(alpha)]

        # Loop through i in alpha getting values for interpolation
        for i in alpha:

            # ------------------------------------------------------
            # Residuals and mu on estimates boundary; dFcHat
            # ------------------------------------------------------

            # Get residuals for field i along dFcHat
            residsi_dFcHat = resids_dFcHat['field'+str(i)]

            # Save residuals for field i along dalphaFc
            resids_dalphaFcHat[str(i)] = residsi_dFcHat[:,dalphaFcHat_loc,:]

            # Get mu for field i along dFcHat
            muHati_dFcHat = muHat_dFcHat['field'+str(i)]

            # Save residuals for field i along dalphaFc
            muHat_dalphaFcHat[str(i)] = muHati_dFcHat[dalphaFcHat_loc,:]

        # Save residuals for alpha
        resids_dFcHat_partitioned[np.array2string(alpha)] = resids_dalphaFcHat

        # Save muhat for alpha
        muHat_dFcHat_partitioned[np.array2string(alpha)] = muHat_dalphaFcHat

    # -------------------------------------------------------------------
    # Get weights from muhat along boundary partitions for interpolation
    # -------------------------------------------------------------------

    # Empty dicts for weights
    weights_dFcHat = {}

    # Loop through boundary partitions
    for alpha in alphas:

        # New empty dicts for weights
        weights_dalphaFcHat = {}

        # Loop through i in alpha getting values for interpolation
        for i in alpha:

            # Get muhati along dalpha Fc
            muHati_dalphaFcHat = muHat_dFcHat_partitioned[np.array2string(alpha)][str(i)]

            # For estimated boundary
            dalphaFcHat_muHati_weights = get_bdry_weights_concat(muHati_dalphaFcHat, c)

            # Save weights
            weights_dalphaFcHat[str(i)] = dalphaFcHat_muHati_weights

        # Save weights
        weights_dFcHat[np.array2string(alpha)] = weights_dalphaFcHat

    # -------------------------------------------------------------------
    # Estimated excursion sets
    # -------------------------------------------------------------------

    # Obtain AcHat for all i
    AcHat = muHats > c

    # Obtain FcHat
    FcHat = cap_muHat > c

    # -------------------------------------------------------------------
    # Perform Bootstrap 
    # -------------------------------------------------------------------
    a_estBdry = bootstrap_resids(resids_dFcHat_partitioned, weights_dFcHat, m, n_boot, p, n_sub)

    # -------------------------------------------------------------------
    # Get FcHat^{+/-}
    # -------------------------------------------------------------------

    # Get the statistic field which defined Achat^{+/-,i}
    g = ((muHats-c)/(sigmas*tau))

    # Take minimum over i
    stat = np.amin(g,axis=0)
    stat = stat.reshape(stat.shape[-2],stat.shape[-1])

    # Obtain FcHat^+ and FcHat^- based on a from the estimated boundary. This variable
    # has axes corresponding to [pvalue, plus/minus, field dimensions]
    FcHat_pm_estBdry = stat >= a_estBdry

    # Return result
    return(FcHat_pm_estBdry[:,0,...], FcHat_pm_estBdry[:,1,...], FcHat, a_estBdry[:,1,...])