import os
import numpy as np
from generateData import *
from boundary import *
from fileio import *

# ===========================================================================
#
# Version 1: Testing for a controlled by
#
# 				a = max(sup_dAc1 |g1|, sup_dAc2 |g2|)
#
# On two cirlces. 
#
# ---------------------------------------------------------------------------
#
# Inputs:
#
# ---------------------------------------------------------------------------
#
# - `OutDir`: Output directory.
# - `nSub`: Number of subjects.
# - `nReals`: Number of realizations.
# - `c`: threshold of interest for mu.
# - `p`: numpy array of p-values.
#
# ---------------------------------------------------------------------------
#
# Developers note: I needed to shorter notation but there wasn't a natural
# 				   shortening for intersection and union so I have used F and
# 				   G for shorthand for intersection and union, respectively,
#				   wherever they occur.
#
# ===========================================================================
def SpatialSims_2mu_v1(OutDir, nSub, nReals, c, p):

    t1overall = time.time()

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

    # Specify signal for each mu
    muSpec1 = {'type': 'circle2D', 'center': np.array([-20,0]), 'fwhm': np.array([5,5]), 'r': 30, 'mag': 3}
    muSpec2 = {'type': 'circle2D', 'center': np.array([20,0]), 'fwhm': np.array([5,5]), 'r': 30, 'mag': 3}

    # Loop through realizations
    for r in np.arange(nReals):

        print('r: ', r)
        # -------------------------------------------------------------------
        # Data generation
        # -------------------------------------------------------------------

        # Obtain data
        data1, mu1 = get_data(muSpec1, data_dim, fwhm)
        data2, mu2 = get_data(muSpec2, data_dim, fwhm)

        # -------------------------------------------------------------------
        # Mean and variance estimates
        # -------------------------------------------------------------------

        # Obtain mu estimate
        muHat1 = np.mean(data1, axis=0).reshape(mu1.shape)
        muHat2 = np.mean(data2, axis=0).reshape(mu1.shape)

        # Obtain sigma
        sigma1 = np.std(data1, axis=0).reshape(mu1.shape)
        sigma2 = np.std(data2, axis=0).reshape(mu1.shape)

		# -------------------------------------------------------------------
        # Boundary locations for Ac1 and Ac2
        # -------------------------------------------------------------------
        # Get boolean maps for the boundary of Ac1 and Ac2
        Ac1_bdry_maps = get_bdry_maps(mu1, c)
        Ac2_bdry_maps = get_bdry_maps(mu2, c)

        # Get coordinates for the boundary of Ac1 and Ac2
        Ac1_bdry_locs = get_bdry_locs(Ac1_bdry_maps)
        Ac2_bdry_locs = get_bdry_locs(Ac2_bdry_maps)

        # Delete maps as we no longer need them
        del Ac1_bdry_maps, Ac2_bdry_maps

        # -------------------------------------------------------------------
        # Boundary locations for Fc (= Ac1 intersect Ac2)
        # -------------------------------------------------------------------
        # Get boolean maps for the boundary of Fc
        Fc_bdry_maps = get_bdry_maps(np.minimum(mu1,mu2), c)

        # Get coordinates for the boundary of Fc
        Fc_bdry_locs = get_bdry_locs(Fc_bdry_maps)

        # Delete maps as we no longer need them
        del Fc_bdry_maps

        # -------------------------------------------------------------------
        # Boundary locations for AcHat1 and AcHat2
        # -------------------------------------------------------------------
        # Get boolean maps for the boundary of AcHat1 and AcHat2
        AcHat1_bdry_maps = get_bdry_maps(muHat1, c)
        AcHat2_bdry_maps = get_bdry_maps(muHat2, c)

        # Get coordinates for the boundary of AcHat1 and AcHat2
        AcHat1_bdry_locs = get_bdry_locs(AcHat1_bdry_maps)
        AcHat2_bdry_locs = get_bdry_locs(AcHat2_bdry_maps)

        # Delete maps as we no longer need them
        del AcHat1_bdry_maps, AcHat2_bdry_maps 

        # -------------------------------------------------------------------
        # Boundary locations for FcHat(= AcHat1 intersect AcHat2)
        # -------------------------------------------------------------------
        # Get boolean maps for the boundary of FcHat
        FcHat_bdry_maps = get_bdry_maps(np.minimum(muHat1,muHat2), c)

        # Get coordinates for the boundary of FcHat
        FcHat_bdry_locs = get_bdry_locs(FcHat_bdry_maps)

        # Delete maps as we no longer need them
        del FcHat_bdry_maps

        # -------------------------------------------------------------------
        # Interpolation weights for Ac1 and Ac2 boundary (Dict version)
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for Ac1 and Ac2
        Ac1_bdry_vals = get_bdry_values(mu1, Ac1_bdry_locs)
        Ac2_bdry_vals = get_bdry_values(mu2, Ac2_bdry_locs)

        # Obtain the weights along the boundary for Ac1 and Ac2
        Ac1_bdry_weights = get_bdry_weights(Ac1_bdry_vals, c)
        Ac2_bdry_weights = get_bdry_weights(Ac2_bdry_vals, c)

        # Delete values as we no longer need them
        del Ac1_bdry_vals, Ac2_bdry_vals

        # -------------------------------------------------------------------
        # Interpolation weights for Fc boundary (Dict version)
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for Fc
        Fc_bdry_vals = get_bdry_values(np.minimum(mu1,mu2), Fc_bdry_locs)

        # Obtain the weights along the boundary for Fc
        Fc_bdry_weights = get_bdry_weights(Fc_bdry_vals, c)

        # Delete values as we no longer need them
        del Fc_bdry_vals

        # -------------------------------------------------------------------
        # Interpolation weights for Ac1 and Ac2 boundary (Array version)
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for Ac1 and Ac2
        Ac1_bdry_vals_concat = get_bdry_values_concat(mu1, Ac1_bdry_locs)
        Ac2_bdry_vals_concat = get_bdry_values_concat(mu2, Ac2_bdry_locs)

        # Obtain the weights along the boundary for Ac1 and Ac2
        Ac1_bdry_weights_concat = get_bdry_weights_concat(Ac1_bdry_vals_concat, c)
        Ac2_bdry_weights_concat = get_bdry_weights_concat(Ac2_bdry_vals_concat, c)

        # Delete values as we no longer need them
        del Ac1_bdry_vals_concat, Ac2_bdry_vals_concat

        # -------------------------------------------------------------------
        # Interpolation weights for Fc boundary (Array version)
        # -------------------------------------------------------------------
		# Obtain the values along the boundary for Fc
        Fc_bdry_vals_concat = get_bdry_values_concat(np.minimum(mu1,mu2), Fc_bdry_locs)

        # Obtain the weights along the boundary for Fc
        Fc_bdry_weights_concat = get_bdry_weights_concat(Fc_bdry_vals_concat, c)
 
 		# Delete values as we no longer need them
        del Fc_bdry_vals_concat
        
        # -------------------------------------------------------------------
        # Interpolation weights for AcHat1 and AcHat2 boundary (Array
        # version)
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for AcHat1 and AcHat2
        AcHat1_bdry_vals_concat = get_bdry_values_concat(muHat1, AcHat1_bdry_locs)
        AcHat2_bdry_vals_concat = get_bdry_values_concat(muHat2, AcHat2_bdry_locs)

        # Obtain the weights along the boundary for AcHat1 and AcHat2
        AcHat1_bdry_weights_concat = get_bdry_weights_concat(AcHat1_bdry_vals_concat, c)
        AcHat2_bdry_weights_concat = get_bdry_weights_concat(AcHat2_bdry_vals_concat, c)

        # Delete values as we no longer need them
        del AcHat1_bdry_vals_concat, AcHat2_bdry_vals_concat

        # -------------------------------------------------------------------
        # Interpolation weights for FcHat boundary (Array version)
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for FcHat
        FcHat_bdry_vals_concat = get_bdry_values_concat(np.minimum(muHat1,muHat2), FcHat_bdry_locs)

        # Obtain the weights along the boundary for FcHat
        FcHat_bdry_weights_concat = get_bdry_weights_concat(FcHat_bdry_vals_concat, c)

        # Delete values as we no longer need them
        del FcHat_bdry_vals_concat

        # -------------------------------------------------------------------
        # Residuals
        # -------------------------------------------------------------------

        # Obtain residuals
        resid1 = (data1-muHat1)/sigma1
        resid2 = (data2-muHat2)/sigma2

		# Residuals along Ac1 and Ac2 boundary
        resid1_Ac1_bdry_concat = get_bdry_values_concat(resid1, Ac1_bdry_locs)
        resid2_Ac2_bdry_concat = get_bdry_values_concat(resid2, Ac2_bdry_locs)

        # Residuals along AcHat1 and AcHat2 boundary
        resid1_AcHat1_bdry_concat = get_bdry_values_concat(resid1, AcHat1_bdry_locs)
        resid2_AcHat2_bdry_concat = get_bdry_values_concat(resid2, AcHat2_bdry_locs)

        # Delete residuals as they are no longer needed
        del resid1, resid2, data1, data2

        # -------------------------------------------------------------------
        # True and estimated excursion sets
        # -------------------------------------------------------------------
        # Obtain Ac1 and Ac2
        Ac1 = mu1 > c
        Ac2 = mu2 > c

        # Obtain Fc
        Fc = np.minimum(mu1,mu2) > c

        # Obtain AcHat2 and AcHat2 
        AcHat1 = muHat1 > c
        AcHat2 = muHat2 > c

        # Obtain FcHat
        FcHat = np.minimum(muHat1,muHat2) > c

        # -------------------------------------------------------------------
        # Muhat1 and MuHat2 (interpolated) along the true Ac1 and Ac2
        # boundary, respectively
        # -------------------------------------------------------------------

        # Residuals along Ac1 and Ac2 boundary
        muHat1_Ac1Bdry = get_bdry_values(muHat1, Ac1_bdry_locs)
        muHat2_Ac2Bdry = get_bdry_values(muHat2, Ac2_bdry_locs)

        # Interpolate along Ac1 and Ac2 boundary
        muHat1_Ac1Bdry = get_bdry_vals_interpolated(muHat1_Ac1Bdry, Ac1_bdry_weights)
        muHat2_Ac2Bdry = get_bdry_vals_interpolated(muHat2_Ac2Bdry, Ac2_bdry_weights)

        # Residuals along Ac1 and Ac2 boundary (array version)
        muHat1_Ac1Bdry_concat = get_bdry_values_concat(muHat1, Ac1_bdry_locs)
        muHat2_Ac2Bdry_concat = get_bdry_values_concat(muHat2, Ac2_bdry_locs)

        # Interpolate along Ac1 and Ac2 boundary (array version)
        muHat1_Ac1Bdry_concat = get_bdry_vals_interpolated_concat(muHat1_Ac1Bdry_concat, Ac1_bdry_weights_concat)
        muHat2_Ac2Bdry_concat = get_bdry_vals_interpolated_concat(muHat2_Ac2Bdry_concat, Ac2_bdry_weights_concat)

        # -------------------------------------------------------------------
        # Bootstrap 
        # -------------------------------------------------------------------
        # Initialize empty bootstrap stores
        max_g1_Ac1 = np.zeros(nBoot)
        max_g2_Ac2 = np.zeros(nBoot)
        max_g1_AcHat1 = np.zeros(nBoot)
        max_g2_AcHat2 = np.zeros(nBoot)

        t1 = time.time()
        # For each bootstrap record the max of the residuals along the
        # boundary
        for b in np.arange(nBoot):

            # Obtain bootstrap variables
            boot_vars = 2*np.random.randint(0,2,boot_dim)-1

            # Reshape for broadcasting purposes (extra axis refers to the fact we have
            # inner and outer boundary values in the last axes of resid_Ac_bdry_concat
            # and resid_AcHat_bdry_concat)
            boot_vars = boot_vars.reshape((*boot_vars.shape),1)

            # Bootstrap residuals along Ac1 and Ac2
            boot_resid1_Ac1_bdry_concat = boot_vars*resid1_Ac1_bdry_concat
            boot_resid2_Ac2_bdry_concat = boot_vars*resid2_Ac2_bdry_concat

            # Bootstrap residuals along AcHat1 and AcHat2
            boot_resid1_AcHat1_bdry_concat = boot_vars*resid1_AcHat1_bdry_concat
            boot_resid2_AcHat2_bdry_concat = boot_vars*resid2_AcHat2_bdry_concat

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of Ac1. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g1_Ac1_bdry_concat = np.zeros(boot_resid1_Ac1_bdry_concat.shape[-2:])
            boot_g1_Ac1_bdry_concat[...,0] = np.sum(boot_resid1_Ac1_bdry_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g1_Ac1_bdry_concat[...,1] = np.sum(boot_resid1_Ac1_bdry_concat[...,1], axis=0)/np.sqrt(nSub)

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of Ac2. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g2_Ac2_bdry_concat = np.zeros(boot_resid2_Ac2_bdry_concat.shape[-2:])
            boot_g2_Ac2_bdry_concat[...,0] = np.sum(boot_resid2_Ac2_bdry_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g2_Ac2_bdry_concat[...,1] = np.sum(boot_resid2_Ac2_bdry_concat[...,1], axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along Ac1. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma1_boot_Ac1_concat = np.zeros(boot_resid1_Ac1_bdry_concat.shape[-2:])
            sigma1_boot_Ac1_concat[...,0] = np.std(boot_resid1_Ac1_bdry_concat[...,0], axis=0, ddof=1)
            sigma1_boot_Ac1_concat[...,1] = np.std(boot_resid1_Ac1_bdry_concat[...,1], axis=0, ddof=1)

            # Obtain bootstrap standard deviations along Ac2. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma2_boot_Ac2_concat = np.zeros(boot_resid2_Ac2_bdry_concat.shape[-2:])
            sigma2_boot_Ac2_concat[...,0] = np.std(boot_resid2_Ac2_bdry_concat[...,0], axis=0, ddof=1)
            sigma2_boot_Ac2_concat[...,1] = np.std(boot_resid2_Ac2_bdry_concat[...,1], axis=0, ddof=1)

            # Divide by the boostrap standard deviation on Ac1
            boot_g1_Ac1_bdry_concat = boot_g1_Ac1_bdry_concat/sigma1_boot_Ac1_concat

            # Divide by the boostrap standard deviation on Ac2
            boot_g2_Ac2_bdry_concat = boot_g2_Ac2_bdry_concat/sigma2_boot_Ac2_concat

			# Sum across subjects to get the bootstrapped a values along
            # the boundary of AcHat1. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g1_AcHat1_bdry_concat = np.zeros(boot_resid1_AcHat1_bdry_concat.shape[-2:])
            boot_g1_AcHat1_bdry_concat[...,0] = np.sum(boot_resid1_AcHat1_bdry_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g1_AcHat1_bdry_concat[...,1] = np.sum(boot_resid1_AcHat1_bdry_concat[...,1], axis=0)/np.sqrt(nSub)

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of AcHat2. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g2_AcHat2_bdry_concat = np.zeros(boot_resid2_AcHat2_bdry_concat.shape[-2:])
            boot_g2_AcHat2_bdry_concat[...,0] = np.sum(boot_resid2_AcHat2_bdry_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g2_AcHat2_bdry_concat[...,1] = np.sum(boot_resid2_AcHat2_bdry_concat[...,1], axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along AcHat1. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma1_boot_AcHat1_concat = np.zeros(boot_resid1_AcHat1_bdry_concat.shape[-2:])
            sigma1_boot_AcHat1_concat[...,0] = np.std(boot_resid1_AcHat1_bdry_concat[...,0], axis=0, ddof=1)
            sigma1_boot_AcHat1_concat[...,1] = np.std(boot_resid1_AcHat1_bdry_concat[...,1], axis=0, ddof=1)

            # Obtain bootstrap standard deviations along AcHat2. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma2_boot_AcHat2_concat = np.zeros(boot_resid2_AcHat2_bdry_concat.shape[-2:])
            sigma2_boot_AcHat2_concat[...,0] = np.std(boot_resid2_AcHat2_bdry_concat[...,0], axis=0, ddof=1)
            sigma2_boot_AcHat2_concat[...,1] = np.std(boot_resid2_AcHat2_bdry_concat[...,1], axis=0, ddof=1)

            # Divide by the boostrap standard deviation on AcHat1
            boot_g1_AcHat1_bdry_concat = boot_g1_AcHat1_bdry_concat/sigma1_boot_AcHat1_concat

            # Divide by the boostrap standard deviation on AcHat2
            boot_g2_AcHat2_bdry_concat = boot_g2_AcHat2_bdry_concat/sigma2_boot_AcHat2_concat

            # Interpolation for Ac1 and Ac2 boundary
            boot_g1_Ac1_bdry_concat = get_bdry_vals_interpolated_concat(boot_g1_Ac1_bdry_concat,Ac1_bdry_weights_concat)
            boot_g2_Ac2_bdry_concat = get_bdry_vals_interpolated_concat(boot_g2_Ac2_bdry_concat,Ac2_bdry_weights_concat)

            # Interpolation for AcHat1 and AcHat2 boundary
            boot_g1_AcHat1_bdry_concat = get_bdry_vals_interpolated_concat(boot_g1_AcHat1_bdry_concat,AcHat1_bdry_weights_concat)
            boot_g2_AcHat2_bdry_concat = get_bdry_vals_interpolated_concat(boot_g2_AcHat2_bdry_concat,AcHat2_bdry_weights_concat)

            # Get maximum along Ac1 and Ac2 boudary
            max_g1_Ac1[b] = np.max(np.abs(boot_g1_Ac1_bdry_concat)) 
            max_g2_Ac2[b] = np.max(np.abs(boot_g2_Ac2_bdry_concat)) 

            # Get maximum along AcHat1 and AcHat2 boudary
            max_g1_AcHat1[b] = np.max(np.abs(boot_g1_AcHat1_bdry_concat))
            max_g2_AcHat2[b] = np.max(np.abs(boot_g2_AcHat2_bdry_concat)) 

        # Get the maximum needed for the true boundary
        max_g_Ac = np.maximum(max_g1_Ac1,max_g2_Ac2)

        # Get the maximum needed for the estimated boundary
        max_g_AcHat = np.maximum(max_g1_AcHat1,max_g2_AcHat2)

        t2 = time.time()
        print('Bootstrap time: ', t2-t1)

        print('max g shape ', max_g_Ac.shape)

        # -------------------------------------------------------------------
        # Obtaining a from percentiles of the max distribution
        # -------------------------------------------------------------------

        # Get the a estimates for the true boundaries and estimated boundaries
        a_trueBdry = np.percentile(max_g_Ac, 100*p).reshape(nPvals,1,1,1)
        a_estBdry = np.percentile(max_g_AcHat, 100*p).reshape(nPvals,1,1,1) # [pvals, 1, [1 for _ in dim>1]]

        # Reformat them to an array form useful for boolean operation
        a_trueBdry = np.concatenate((-a_trueBdry,a_trueBdry),axis=1)
        a_estBdry = np.concatenate((-a_estBdry,a_estBdry),axis=1)

        # Get the statistic field which defined Achat^{+/-,1/2}
        stat1 = ((muHat1-c)/(sigma1*tau)).reshape(1,(*muHat1.shape))
        stat2 = ((muHat2-c)/(sigma2*tau)).reshape(1,(*muHat2.shape))
        
        # Minimum for intersection
        stat = np.minimum(stat1,stat2)
        stat = stat.reshape(stat.shape[-2],stat.shape[-1])

        # Obtain FcHat^+ and FcHat^- based on a from the true boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        FcHat_pm_trueBdry = stat >= a_trueBdry

        # Obtain FcHat^+ and FcHat^- based on a from the estimated boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        FcHat_pm_estBdry = stat >= a_estBdry

        # -------------------------------------------------------------------
        # Some set logic to work out violations
        # -------------------------------------------------------------------

        # Obtain FcHat^+\FcHat based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        FcHatp_sub_Fc_trueBdry = FcHat_pm_trueBdry[:,1,...] & ~Fc[...]

        # Obtain FcHat^+\FcHat based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        FcHatp_sub_Fc_estBdry = FcHat_pm_estBdry[:,1,...] & ~Fc[...]

        # Obtain Fc\FcHat^- based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        Fc_sub_FcHatm_trueBdry = Fc[...] & ~FcHat_pm_trueBdry[:,0,...]

        # Obtain Fc\FcHat^- based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        Fc_sub_FcHatm_estBdry = Fc[...] & ~FcHat_pm_estBdry[:,0,...]

        print(Fc_sub_FcHatm_estBdry.shape, Fc.shape, FcHat_pm_estBdry.shape)

        # -------------------------------------------------------------------
        # Check whether there were any boundary violations using voxelwise
        # set logic (checking if voxels existed in one set but not another, 
        # etc)
        # -------------------------------------------------------------------
        # Record if we saw a violation in the true boundary based sets
        trueBdry_success[r,:] = 1-(np.any(FcHatp_sub_Fc_trueBdry,axis=(1,2)) | np.any(Fc_sub_FcHatm_trueBdry,axis=(1,2))) # MARKER: AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success[r,:] = 1-(np.any(FcHatp_sub_Fc_estBdry,axis=(1,2)) | np.any(Fc_sub_FcHatm_estBdry,axis=(1,2)))# MARKER: AXES WONT WORK FOR 3D ATM

        # -------------------------------------------------------------------
        # Get stat along the Fc boundary
        # -------------------------------------------------------------------
        # Get the values along the outer and inner boundaries
        stat_FcBdry = get_bdry_values(stat, Fc_bdry_locs)

        # Interpolate to get the values along the true boundary
        stat_FcBdry = get_bdry_vals_interpolated(stat_FcBdry, Fc_bdry_weights)

        # -------------------------------------------------------------------
        # Check whether there were any boundary violations using interpolated
        # boundary values (checking if voxels had values corresponding to no
        # violations, etc)
        # -------------------------------------------------------------------

        # Perform lower check on stat map using thresholds based on the
        # estimated boundary
        bdry_lowerCheck_estBdry = stat_FcBdry >= a_estBdry[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # estimated boundary
        bdry_upperCheck_estBdry = stat_FcBdry <= a_estBdry[:,1,:,0]

        # Perform lower check on stat map using thresholds based on the
        # true boundary
        bdry_lowerCheck_trueBdry = stat_FcBdry >= a_trueBdry[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # true boundary
        bdry_upperCheck_trueBdry = stat_FcBdry <= a_trueBdry[:,1,:,0]

        print('here')
        print(bdry_upperCheck_trueBdry.shape)
        print(stat_FcBdry.shape)
        print(a_trueBdry[:,1,:,0].shape)

        # -------------------------------------------------------------------
        # Work out whether simulation observed successful sets.
        # -------------------------------------------------------------------
        # Record if we saw a violation in the true boundary based sets
        trueBdry_success_intrp[r,:] = np.all(bdry_lowerCheck_trueBdry,axis=(1)) & np.all(bdry_upperCheck_trueBdry,axis=(1))# MARKER: AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success_intrp[r,:] = np.all(bdry_lowerCheck_estBdry,axis=(1)) & np.all(bdry_upperCheck_estBdry,axis=(1)) # MARKER: AXES WONT WORK FOR 3D ATM


    # For the interpolated boundary success checks, we still need to do the 
    # voxelwise checks as well. This will take care of that.
    trueBdry_success_intrp = trueBdry_success_intrp*trueBdry_success
    estBdry_success_intrp = estBdry_success_intrp*estBdry_success

    # Coverage probabilities
    coverage_trueBdry = np.mean(trueBdry_success,axis=0)
    coverage_estBdry = np.mean(estBdry_success,axis=0)

    # Coverage probabilities
    coverage_trueBdry_intrp = np.mean(trueBdry_success_intrp,axis=0)
    coverage_estBdry_intrp = np.mean(estBdry_success_intrp,axis=0)

    print('Coverage: ', coverage_estBdry_intrp)

    # Save the violations to a file
    append_to_file('trueSuccess'+str(nSub)+'v1.csv', trueBdry_success) 
    append_to_file('estSuccess'+str(nSub)+'v1.csv', estBdry_success)
    append_to_file('trueSuccess'+str(nSub)+'v1_intrp.csv', trueBdry_success_intrp) 
    append_to_file('estSuccess'+str(nSub)+'v1_intrp.csv', estBdry_success_intrp)

    t2overall = time.time()

    print('overall time: ', t2overall-t1overall)

SpatialSims_2mu_v1('/home/tommaullin/Documents/ConfSets/',100, 30, 2, np.linspace(0,1,21))