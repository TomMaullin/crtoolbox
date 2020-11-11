import os
import numpy as np
from generateData import *
from boundary import *
from fileio import *

# ===========================================================================
#
# Inputs:
#
# ---------------------------------------------------------------------------
#
# - `OutDir`: Output directory.
# - `nSub`: Number of subjects.
# - `simType`: String representing simulation type. e.g. `ramp` for ramp.
# - `nReals`: Number of realizations.
# - `c`: threshold of interest for mu.
#
# ===========================================================================
def SpatialSims(OutDir, nSub, simType, nReals, c, p):

    # Define tau_n
    tau = 1/np.sqrt(nSub)

    # Get the number of p-values we're looking at
    nPvals = len(p)

    # Define the number of bootstraps
    nBoot = 5000 # Recommended 1e4

    # Dimensions of simulated data
    data_dim = np.array([nSub, 100,100])

    # Dimensions of bootstrap variables
    boot_dim = np.array([nSub, 1]) # MARKER: COULD MAKE  np.array([batchBoot, nSub, 1])

    # Smoothing
    fwhm = [0,3,3]

    # Initialise array for recording the whether set violations occured, for a derived
    # from bootstrapping the true boundary and estimated boundary, respectively. The
    # result in these arrays are based on voxelwise assessment of set condition 
    # violations.
    trueBdry_success = np.zeros((nReals,nPvals))
    estBdry_success = np.zeros((nReals,nPvals))

    # Initialise array for recording the whether set violations occured, for a derived
    # from bootstrapping the true boundary and estimated boundary, respectively. The
    # result in these arrays are based on interpolation based assessment of set
    # condition violations.
    trueBdry_success_intrp = np.zeros((nReals,nPvals))
    estBdry_success_intrp = np.zeros((nReals,nPvals))

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
        # Muhat (interpolated) along the true Ac boundary
        # -------------------------------------------------------------------

        # Residuals along Ac boundary
        muHat_AcBdry = get_bdry_values(muHat, Ac_bdry_locs)

        # Interpolate along Ac boundary
        muHat_AcBdry = get_bdry_vals_interpolated(muHat_AcBdry, Ac_bdry_weights)

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
        a_trueBdry = np.percentile(max_g_Ac, 100*p).reshape(nPvals,1,1,1)
        a_estBdry = np.percentile(max_g_AcHat, 100*p).reshape(nPvals,1,1,1) # [pvals, 1, [1 for _ in dim>1]]

        # Reformat them to an array form useful for boolean operation
        a_trueBdry = np.concatenate((-a_trueBdry,a_trueBdry),axis=1)
        a_estBdry = np.concatenate((-a_estBdry,a_estBdry),axis=1)

        # Get the statistic field which defined Achat^{+/-}
        stat = ((muHat-c)/(sigma*tau)).reshape(1,(*muHat.shape))

        # Obtain AcHat^+ and AcHat^- based on a from the true boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        AcHat_pm_trueBdry = stat >= a_trueBdry

        # Obtain AcHat^+ and AcHat^- based on a from the estimated boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        AcHat_pm_estBdry = stat >= a_estBdry

        # -------------------------------------------------------------------
        # Some set logic to work out violations
        # -------------------------------------------------------------------

        # Obtain AcHat^+\AcHat based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        AcHatp_sub_Ac_trueBdry = AcHat_pm_trueBdry[:,1,...] & ~Ac[...]

        # Obtain AcHat^+\AcHat based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        AcHatp_sub_Ac_estBdry = AcHat_pm_estBdry[:,1,...] & ~Ac[...]

        # Obtain Ac\AcHat^- based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        Ac_sub_AcHatm_trueBdry = Ac[...] & ~AcHat_pm_trueBdry[:,0,...]

        # Obtain Ac\AcHat^- based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        Ac_sub_AcHatm_estBdry = Ac[...] & ~AcHat_pm_estBdry[:,0,...]

        # -------------------------------------------------------------------
        # Check whether there were any boundary violations using voxelwise
        # set logic (checking if voxels existed in one set but not another, 
        # etc)
        # -------------------------------------------------------------------

        print('here')
        print(Ac_sub_AcHatm_estBdry.shape)
        print(Ac_sub_AcHatm_trueBdry.shape)

        # Record if we saw a violation in the true boundary based sets
        trueBdry_success[r,:] = 1-(np.any(AcHatp_sub_Ac_trueBdry,axis=(1,2)) | np.any(Ac_sub_AcHatm_trueBdry,axis=(1,2))) # MARKER: AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success[r,:] = 1-(np.any(AcHatp_sub_Ac_estBdry,axis=(1,2)) | np.any(Ac_sub_AcHatm_estBdry,axis=(1,2)))# MARKER: AXES WONT WORK FOR 3D ATM

        # -------------------------------------------------------------------
        # Work out the threshold values of the muhat image given by:
        #               thr(s) = c +/- a\tau\sigma(s)
        # -------------------------------------------------------------------

        # Obtain threshold maps for muhat to get lower and upper sets (for 
        # true boundary). This variable will have dimensions corresponding
        # to [number of p values, -a or a, dimensions of mu].
        muHat_threshs_trueBdry = c + a_trueBdry*tau*sigma

        # Obtain threshold maps for muhat to get lower and upper sets (for 
        # estimated boundary). This variable will have dimensions
        # corresponding to [number of p values, -a or a, dimensions of mu].
        muHat_threshs_estBdry = c + a_estBdry*tau*sigma

        # -------------------------------------------------------------------
        # Work out the thresholds along the true boundary, using `a` derived
        # from the maxima along the true boundary,
        # -------------------------------------------------------------------
        # Thresholds for muhat along Ac boundary
        muHat_AcBdry_threshs_trueBdry = get_bdry_values(muHat_threshs_trueBdry, Ac_bdry_locs)

        # Delete whole map as we no longer need it
        del muHat_threshs_trueBdry

        # Interpolate along Ac boundary
        muHat_AcBdry_threshs_trueBdry = get_bdry_vals_interpolated(muHat_AcBdry_threshs_trueBdry, Ac_bdry_weights)

        # -------------------------------------------------------------------
        # Work out the thresholds along the true boundary, using `a` derived
        # from the maxima along the estimated boundary,
        # -------------------------------------------------------------------
        # Thresholds for muhat along Ac boundary
        muHat_AcBdry_threshs_estBdry = get_bdry_values(muHat_threshs_estBdry, Ac_bdry_locs)

        # Delete whole map as we no longer need it
        del muHat_threshs_estBdry

        # Interpolate along Ac boundary
        muHat_AcBdry_threshs_estBdry = get_bdry_vals_interpolated(muHat_AcBdry_threshs_estBdry, Ac_bdry_weights)

        print(muHat_AcBdry_threshs_estBdry.shape, muHat_AcBdry_threshs_trueBdry.shape)

        # -------------------------------------------------------------------
        # Check whether there were any boundary violations using interpolated
        # boundary values (checking if voxels had values corresponding to no
        # violations, etc)
        # -------------------------------------------------------------------

        print('marker')
        print(muHat_AcBdry_threshs_trueBdry)

        print(muHat_AcBdry_threshs_trueBdry.shape)

        # Perform lower check on muHat using thresholds based on the
        # estimated boundary
        bdry_lowerCheck_estBdry = muHat_AcBdry >= muHat_AcBdry_threshs_estBdry[:,0,...]

        # Perform upper check on muHat using thresholds based on the
        # estimated boundary
        bdry_upperCheck_estBdry = muHat_AcBdry <= muHat_AcBdry_threshs_estBdry[:,1,...]

        # Perform lower check on muHat using thresholds based on the
        # true boundary
        bdry_lowerCheck_trueBdry = muHat_AcBdry >= muHat_AcBdry_threshs_trueBdry[:,0,...]

        # Perform upper check on muHat using thresholds based on the
        # true boundary
        bdry_upperCheck_trueBdry = muHat_AcBdry <= muHat_AcBdry_threshs_trueBdry[:,1,...]

        print('check shape')
        print(bdry_lowerCheck_trueBdry.shape)

        # -------------------------------------------------------------------
        # Work out whether simulation observed successful sets.
        # -------------------------------------------------------------------
        # Record if we saw a violation in the true boundary based sets
        trueBdry_success_intrp[r,:] = np.all(bdry_lowerCheck_trueBdry,axis=(1)) & np.all(bdry_upperCheck_trueBdry,axis=(1))# MARKER: AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success_intrp[r,:] = np.all(bdry_lowerCheck_estBdry,axis=(1)) & np.all(bdry_upperCheck_estBdry,axis=(1)) # MARKER: AXES WONT WORK FOR 3D ATM

    # For the interpolated boundary success checks, we still need to do the 
    # voxelwise checks as well. This will take care of that.
    print('successes')
    print(trueBdry_success_intrp.shape, trueBdry_success.shape)
    print(trueBdry_success_intrp, trueBdry_success)
    trueBdry_success_intrp = trueBdry_success_intrp*trueBdry_success
    estBdry_success_intrp = estBdry_success_intrp*estBdry_success

    # Coverage probabilities
    coverage_trueBdry = np.mean(trueBdry_success,axis=0)
    coverage_estBdry = np.mean(estBdry_success,axis=0)

    # Coverage probabilities
    coverage_trueBdry_intrp = np.mean(trueBdry_success_intrp,axis=0)
    coverage_estBdry_intrp = np.mean(estBdry_success_intrp,axis=0)

    print(coverage_estBdry_intrp, coverage_estBdry)

    # Save the violations to a file
    append_to_file('trueSuccess'+str(nSub)+'.csv', trueBdry_success) 
    append_to_file('estSuccess'+str(nSub)+'.csv', estBdry_success)
    append_to_file('trueSuccess'+str(nSub)+'_intrp.csv', trueBdry_success_intrp) 
    append_to_file('estSuccess'+str(nSub)+'_intrp.csv', estBdry_success_intrp)

# Run example
SpatialSims('/home/tommaullin/Documents/ConfSets/',100, 'rampHoriz2D', 30, 2, np.linspace(0,1,21))