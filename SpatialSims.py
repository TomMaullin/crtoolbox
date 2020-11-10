import os
import numpy as np
from generateData import *
from boundary import *

# ===========================================================================
#
# Inputs:
#
# ---------------------------------------------------------------------------
#
# - `nSub`: Number of subjects.
# - `simType`: String representing simulation type. e.g. `ramp` for ramp.
# - `nReals`: Number of realizations.
# - `c`: threshold of interest for mu.
#
# ===========================================================================
def SpatialSims(nSub, simType, nReals, c, p):

    # Define tau_n
    tau = 1/np.sqrt(nSub)

    # Define the number of bootstraps
    nBoot = 5000 # Recommended 1e4

    # Dimensions of simulated data
    data_dim = np.array([nSub, 100,100])

    # Dimensions of bootstrap variables
    boot_dim = np.array([nSub, 1]) # MARKER: COULD MAKE  np.array([batchBoot, nSub, 1])

    # Smoothing
    fwhm = [0,3,3]

    # Loop through realizations
    for r in np.arange(nReals):

        print('r: ', r)
        # -------------------------------------------------------------------
        # Data generation
        # -------------------------------------------------------------------

        # Obtain data
        data, mu = get_data(simType, data_dim, fwhm)

        # -------------------------------------------------------------------
        # Mean and variance estimates
        # -------------------------------------------------------------------

        # Obtain mu estimate
        muHat = np.mean(data, axis=0).reshape(mu.shape)

        # Obtain sigma
        sigma = np.std(data, axis=0).reshape(mu.shape)

        # -------------------------------------------------------------------
        # Boundary locations for Ac
        # -------------------------------------------------------------------
        # Get boolean maps for the boundary of Ac
        Ac_bdry_maps = get_bdry_maps(mu, c)

        # Get coordinates for the boundary of Ac
        Ac_bdry_locs = get_bdry_locs(Ac_bdry_maps)

        # Delete maps as we no longer need them
        del Ac_bdry_maps

        # -------------------------------------------------------------------
        # Boundary locations for AcHat
        # -------------------------------------------------------------------
        # Get boolean maps for the boundary of AcHat
        AcHat_bdry_maps = get_bdry_maps(muHat, c)

        # Get coordinates for the boundary of AcHat
        AcHat_bdry_locs = get_bdry_locs(AcHat_bdry_maps)

        # Delete maps as we no longer need them
        del AcHat_bdry_maps

        # -------------------------------------------------------------------
        # Interpolation weights for Ac boundary
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for Ac
        Ac_bdry_vals = get_bdry_values(mu, Ac_bdry_locs)

        # Obtain the weights along the boundary for Ac
        Ac_bdry_weights = get_bdry_weights(Ac_bdry_vals, c)

        # Delete values as we no longer need them
        del Ac_bdry_vals

        # -------------------------------------------------------------------
        # Interpolation weights for AcHat boundary
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for AcHat
        AcHat_bdry_vals = get_bdry_values(muHat, AcHat_bdry_locs)

        # Obtain the weights along the boundary for Ac
        AcHat_bdry_weights = get_bdry_weights(AcHat_bdry_vals, c)

        # Delete values as we no longer need them
        del AcHat_bdry_vals

        # -------------------------------------------------------------------
        # Residuals
        # -------------------------------------------------------------------

        # Obtain residuals
        resid = (muHat-data)/sigma

        # Residuals along Ac boundary
        resid_Ac_bdry = get_bdry_values(resid, Ac_bdry_locs)

        # Interpolate along Ac boundary
        resid_Ac_bdry = get_bdry_vals_interpolated(resid_Ac_bdry, Ac_bdry_weights)

        # Residuals along AcHat boundary
        resid_AcHat_bdry = get_bdry_values(resid, AcHat_bdry_locs)

        # Interpolate along AcHat boundary
        resid_AcHat_bdry = get_bdry_vals_interpolated(resid_AcHat_bdry, AcHat_bdry_weights)

        # Delete residuals as they are no longer needed
        del resid

        # -------------------------------------------------------------------
        # True and estimated excursion sets
        # -------------------------------------------------------------------
        # Obtain Ac
        Ac = mu > c

        # Obtain estimated Ac
        AcHat = muHat > c

        # -------------------------------------------------------------------
        # Bootstrap # MARKER: CAN DEFO DELETE THIS LOOP
        # -------------------------------------------------------------------

        # Initialize empty bootstrap stores
        max_g_Ac = np.zeros(nBoot)
        max_g_AcHat = np.zeros(nBoot)

        t1 = time.time()
        # For each bootstrap record the max of the residuals along the
        # boundary
        for b in np.arange(nBoot):

            #print('b: ', b)

            # Obtain bootstrap variables
            boot_vars = 2*np.random.randint(0,2,boot_dim)-1

            # Bootstrap residuals along Ac
            boot_resid_Ac_bdry = boot_vars*resid_Ac_bdry

            # Bootstrap residuals along AcHat
            boot_resid_AcHat_bdry = boot_vars*resid_AcHat_bdry

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of Ac
            boot_g_Ac_bdry = np.sum(boot_resid_Ac_bdry, axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along Ac
            sigma_boot_Ac = np.std(boot_resid_Ac_bdry, axis=0)

            # Divide by the boostrap standard deviation on Ac
            boot_g_Ac_bdry = boot_g_Ac_bdry/sigma_boot_Ac
            
            # Sum across subjects to get the bootstrapped g values along
            # the boundary of AcHat
            boot_g_AcHat_bdry = np.sum(boot_resid_AcHat_bdry, axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along AcHat
            sigma_boot_AcHat = np.std(boot_resid_AcHat_bdry, axis=0)

            # Divide by the boostrap standard deviation on AcHat
            boot_g_AcHat_bdry = boot_g_AcHat_bdry/sigma_boot_AcHat

            # Get maximum along Ac boudary
            max_g_Ac[b] = np.max(boot_g_Ac_bdry) 

            # Get maximum along AcHat boudary
            max_g_AcHat[b] = np.max(boot_g_AcHat_bdry) 

        t2 = time.time()
        print('Bootstrap time: ', t2-t1)
        print(max_g_Ac)
        print(max_g_AcHat)

        # -------------------------------------------------------------------
        # Obtaining a from percentiles of the max distribution
        # -------------------------------------------------------------------

        # Get the a estimates along the true Ac boundary and estimated AcHat boundary
        a_trueBdry = np.percentile(max_g_Ac, 100*p)
        a_estBdry = np.percentile(max_g_AcHat, 100*p) # MARKER: CAN MAKE THIS WORK FOR MULTI P - JUST NEED TO LOOK INTO IT

        # Reformat them to an array form useful for boolean operation
        a_trueBdry = np.array([-a_trueBdry,a_trueBdry]).reshape(2,1,1)
        a_estBdry = np.array([-a_estBdry,a_estBdry]).reshape(2,1,1)

        # Obtain AcHat^+ and AcHat^- based on a from the true boundary
        Ac_pm_trueBdry = ((muHat-c)/(sigma*tau)) >= a_trueBdry

        # Obtain AcHat^+ and AcHat^- based on a from the estimated boundary
        Ac_pm_estBdry = ((muHat-c)/(sigma*tau)) >= a_estBdry

SpatialSims(100, 'rampHoriz2D', 5000, 2, 0.95)