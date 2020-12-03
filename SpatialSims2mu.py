import os
import numpy as np
from generateData import *
from boundary import *
from setTheory2mu import *
from fileio import *

# ===========================================================================
#
# Inputs:
#
# ---------------------------------------------------------------------------
#
# - `OutDir`: Output directory.
# - `nSub`: Number of subjects.
# - `muSpec`: Dictionary specifying mu to simulate. Always must include a 
#             `type` parameter specifying `ramp2D` or `circle2D`.
#             ----------------------------------------------------------------
#             For a ramp the following must also be given:
#                - a: lower value of ramp
#                - b: upper value of ramp
#                - orient: string representing `vertical` or `horizontal`
#             ----------------------------------------------------------------
#             For a circle the following must also be given:
#                - center: center coordinates for circle (treating origin
#                          as center of grid)
#                - r: radius of circle
#                - fwhm: fwhm used to smooth the circle 
#                - mag: Magnitude of signal
#             ----------------------------------------------------------------
# - `nReals`: Number of realizations.
# - `c`: threshold of interest for mu.
# - `p`: numpy array of p-values.
#
# ===========================================================================
def SpatialSims2mu(OutDir, nSub, muSpec, nReals, c, p):

    t1overall = time.time()

    # Number of mus
    nMu = 2

    # Define tau_n
    tau = 1/np.sqrt(nSub)

    # Get the number of p-values we're looking at
    nPvals = len(p)

    # Define the number of bootstraps
    nBoot = 5000 # Recommended 1e4

    # Dimensions of simulated data
    data_dim = np.array([nSub,nMu,100,100])

    # Dimensions of bootstrap variables
    boot_dim = np.array([nSub, 1]) # MARKER: COULD MAKE  np.array([batchBoot, nSub, 1])

    # Smoothing
    fwhm = [0,0,3,3]

    # Initialise array for recording the whether set violations occured, for a derived
    # from bootstrapping the true boundary and estimated boundary, respectively. The
    # result in these arrays are based on voxelwise assessment of set condition 
    # violations.
    trueBdry_success = np.zeros((nReals,nPvals))
    estBdry_success = np.zeros((nReals,nPvals))

    # Repeat for union boundary
    trueBdry_success_cup = np.zeros((nReals,nPvals))
    estBdry_success_cup = np.zeros((nReals,nPvals))

    # Repeat for intersection boundary
    trueBdry_success_cap = np.zeros((nReals,nPvals))
    estBdry_success_cap = np.zeros((nReals,nPvals))

    # Initialise array for recording the whether set violations occured, for a derived
    # from bootstrapping the true boundary and estimated boundary, respectively. The
    # result in these arrays are based on interpolation based assessment of set
    # condition violations.
    trueBdry_success_intrp = np.zeros((nReals,nPvals))
    estBdry_success_intrp = np.zeros((nReals,nPvals))

    # Repeat for union boundary
    trueBdry_success_intrp_cup = np.zeros((nReals,nPvals))
    estBdry_success_intrp_cup = np.zeros((nReals,nPvals))

    # Repeat for intersection boundary
    trueBdry_success_intrp_cap = np.zeros((nReals,nPvals))
    estBdry_success_intrp_cap = np.zeros((nReals,nPvals))

    # Loop through realizations
    for r in np.arange(nReals):

        print('r: ', r)
        # -------------------------------------------------------------------
        # Data generation
        # -------------------------------------------------------------------

        # Obtain data
        data, mu = get_data(muSpec, data_dim, fwhm)

        # Reshape slightly
        mu = mu.reshape(nMu,100,100)

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

        # Get the boundary locations
        Ac1_bdry_locs, Ac2_bdry_locs, Ac_cap_bdry_locs, Ac_cup_bdry_locs, Ac_bdry_cap_locs, Ac_bdry_cup_locs = get_2mu_bdry_locs(mu,c,image=False)

        # -------------------------------------------------------------------
        # Boundary locations for AcHat
        # -------------------------------------------------------------------

        # Get the boundary locations
        AcHat1_bdry_locs, AcHat2_bdry_locs, AcHat_cap_bdry_locs, AcHat_cup_bdry_locs, AcHat_bdry_cap_locs, AcHat_bdry_cup_locs = get_2mu_bdry_locs(muHat,c,image=False)

        print('marker')
        # -------------------------------------------------------------------
        # Interpolation weights for Ac1 boundary
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for Ac1
        Ac1_bdry_vals = get_bdry_values(mu, Ac1_bdry_locs)

        # Obtain the weights along the boundary for Ac1
        Ac1_bdry_weights = get_bdry_weights(Ac1_bdry_vals, c)

        # Delete values as we no longer need them
        del Ac1_bdry_vals

        # -------------------------------------------------------------------
        # Interpolation weights for Ac2 boundary
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for Ac2
        Ac2_bdry_vals = get_bdry_values(mu, Ac2_bdry_locs)

        # Obtain the weights along the boundary for Ac2
        Ac2_bdry_weights = get_bdry_weights(Ac2_bdry_vals, c)

        # Delete values as we no longer need them
        del Ac2_bdry_vals

        # -------------------------------------------------------------------
        # Interpolation weights for union boundary
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for Ac union
        Ac_cup_bdry_vals = get_bdry_values(mu, Ac_cup_bdry_locs)

        # Obtain the weights along the boundary for Ac union
        Ac_cup_bdry_weights = get_bdry_weights(Ac_cup_bdry_vals, c)

        # Delete values as we no longer need them
        del Ac_cup_bdry_vals

        # -------------------------------------------------------------------
        # Interpolation Ac1, Ac2 weights for intersection boundary
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for Ac intersection
        Ac_cap_bdry_vals = get_bdry_values(mu, Ac_cap_bdry_locs)

        # Obtain the weights along the boundary for Ac intersection
        Ac_cap_bdry_weights = get_bdry_weights(Ac_cap_bdry_vals, c)

        # Delete values as we no longer need them
        del Ac_cap_bdry_vals

        # -------------------------------------------------------------------
        # Interpolation Fc weights for intersection boundary 
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for Ac intersection
        Fc_bdry_vals = get_bdry_values(np.amin(mu,axis=0), Ac_cap_bdry_locs)

        # Obtain the weights along the boundary for Ac intersection
        Fc_bdry_weights = get_bdry_weights(Fc_bdry_vals, c)

        # Delete values as we no longer need them
        del Fc_bdry_vals

        # -------------------------------------------------------------------
        # Get intersection and union fields for mu and muhat 
        # -------------------------------------------------------------------
        mu_cup, mu_cap = get_2mu_cmbtns(mu)
        muHat_cup, muHat_cap = get_2mu_cmbtns(muHat)

        # ---------------------------------------------------------------
        # Interpolation weights for Ac1 boundary
        # ---------------------------------------------------------------
        # Obtain the values along the boundary for Ac1
        Ac1_bdry_vals_concat = get_bdry_values_concat(mu[0,:,:], Ac1_bdry_locs)

        # Obtain the weights along the boundary for Ac1
        Ac1_bdry_weights_concat = get_bdry_weights_concat(Ac1_bdry_vals_concat, c)

        # Delete values as we no longer need them
        del Ac1_bdry_vals_concat

        # ---------------------------------------------------------------
        # Interpolation weights for Ac2 boundary
        # ---------------------------------------------------------------
        # Obtain the values along the boundary for Ac2
        Ac2_bdry_vals_concat = get_bdry_values_concat(mu[1,:,:], Ac2_bdry_locs)

        # Obtain the weights along the boundary for Ac2
        Ac2_bdry_weights_concat = get_bdry_weights_concat(Ac2_bdry_vals_concat, c)

        # Delete values as we no longer need them
        del Ac2_bdry_vals_concat

        # ---------------------------------------------------------------
        # Interpolation weights for union Ac boundary
        # ---------------------------------------------------------------
        # Obtain the values along the boundary for Ac union
        Ac_cup_bdry_vals_concat = get_bdry_values_concat(mu_cup, Ac_cup_bdry_locs)

        # Obtain the weights along the boundary for Ac union
        Ac_cup_bdry_weights_concat = get_bdry_weights_concat(Ac_cup_bdry_vals_concat, c)

        # Delete values as we no longer need them
        del Ac_cup_bdry_vals_concat

        # ---------------------------------------------------------------
        # Interpolation weights for intersection Ac boundary
        # ---------------------------------------------------------------
        # Obtain the values along the boundary for Ac intersection
        Ac_cap_bdry_vals_concat = get_bdry_values_concat(mu_cap, Ac_cap_bdry_locs)

        # Obtain the weights along the boundary for Ac intersection
        Ac_cap_bdry_weights_concat = get_bdry_weights_concat(Ac_cap_bdry_vals_concat, c)

        # Delete values as we no longer need them
        del Ac_cap_bdry_vals_concat

        # ---------------------------------------------------------------
        # Interpolation weights for AcHat1 boundary
        # ---------------------------------------------------------------
        # Obtain the values along the boundary for AcHat1
        AcHat1_bdry_vals_concat = get_bdry_values_concat(mu[0,:,:], AcHat1_bdry_locs)

        print('feiwfiwnv', AcHat1_bdry_vals_concat.shape)
        print(np.max(AcHat1_bdry_vals_concat))
        print(np.min(AcHat1_bdry_vals_concat))

        # Obtain the weights along the boundary for AcHat1
        AcHat1_bdry_weights_concat = get_bdry_weights_concat(AcHat1_bdry_vals_concat, c)

        print(AcHat1_bdry_weights_concat.shape)
        print(np.max(AcHat1_bdry_weights_concat))
        # Delete values as we no longer need them
        del AcHat1_bdry_vals_concat

        # ---------------------------------------------------------------
        # Interpolation weights for AcHat2 boundary
        # ---------------------------------------------------------------
        # Obtain the values along the boundary for AcHat2
        AcHat2_bdry_vals_concat = get_bdry_values_concat(mu[1,:,:], AcHat2_bdry_locs)

        # Obtain the weights along the boundary for AcHat2
        AcHat2_bdry_weights_concat = get_bdry_weights_concat(AcHat2_bdry_vals_concat, c)

        # Delete values as we no longer need them
        del AcHat2_bdry_vals_concat

        # ---------------------------------------------------------------
        # Interpolation weights for union AcHat boundary
        # ---------------------------------------------------------------
        # Obtain the values along the boundary for AcHat union
        AcHat_cup_bdry_vals_concat = get_bdry_values_concat(mu_cup, AcHat_cup_bdry_locs)

        # Obtain the weights along the boundary for AcHat union
        AcHat_cup_bdry_weights_concat = get_bdry_weights_concat(AcHat_cup_bdry_vals_concat, c)

        # Delete values as we no longer need them
        del AcHat_cup_bdry_vals_concat

        # ---------------------------------------------------------------
        # Interpolation weights for intersection AcHat boundary
        # ---------------------------------------------------------------
        # Obtain the values along the boundary for AcHat intersection
        AcHat_cap_bdry_vals_concat = get_bdry_values_concat(mu_cap, AcHat_cap_bdry_locs)

        # Obtain the weights along the boundary for AcHat intersection
        AcHat_cap_bdry_weights_concat = get_bdry_weights_concat(AcHat_cap_bdry_vals_concat, c)

        # Delete values as we no longer need them
        del AcHat_cap_bdry_vals_concat

        # -------------------------------------------------------------------
        # Residuals
        # -------------------------------------------------------------------

        # Obtain residuals
        resid = (data-muHat)/sigma

        # Residuals for mu1 along Ac1 boundary
        resid1_Ac1_bdry_concat = get_bdry_values_concat(resid[:,0,:,:], Ac1_bdry_locs)

        # Residuals for mu2 along Ac2 boundary
        resid2_Ac2_bdry_concat = get_bdry_values_concat(resid[:,1,:,:], Ac2_bdry_locs)

        # Residuals for both mu along Ac union boundary
        resid_Ac_cup_bdry_concat = get_bdry_values_concat(resid[:,:,:,:], Ac_cup_bdry_locs)

        # Residuals for both mu along Ac intersection boundary
        resid_Ac_cap_bdry_concat = get_bdry_values_concat(resid[:,:,:,:], Ac_cap_bdry_locs)

        # Residuals for mu1 along AcHat1 boundary
        resid1_AcHat1_bdry_concat = get_bdry_values_concat(resid[:,0,:,:], AcHat1_bdry_locs)

        # Residuals for mu2 along AcHat2 boundary
        resid2_AcHat2_bdry_concat = get_bdry_values_concat(resid[:,1,:,:], AcHat2_bdry_locs)

        # Residuals for both mu along AcHat union boundary
        resid_AcHat_cup_bdry_concat = get_bdry_values_concat(resid[:,:,:,:], AcHat_cup_bdry_locs)

        # Residuals for both mu along AcHat intersection boundary
        resid_AcHat_cap_bdry_concat = get_bdry_values_concat(resid[:,:,:,:], AcHat_cap_bdry_locs)

        # print('shapes: ')
        # print(resid1_Ac1_bdry_concat.shape)
        # print(resid2_Ac2_bdry_concat.shape)
        # print(resid_cup_Ac_cup_bdry_concat.shape)
        # print(resid_cap_Ac_cap_bdry_concat.shape)
        # print(resid1_AcHat1_bdry_concat.shape)
        # print(resid2_AcHat2_bdry_concat.shape)
        # print(resid_cup_AcHat_cup_bdry_concat.shape)
        # print(resid_cap_AcHat_cap_bdry_concat.shape)


        # Delete residuals as they are no longer needed
        del resid, data

        # -------------------------------------------------------------------
        # True and estimated excursion sets
        # -------------------------------------------------------------------
        # Obtain Ac
        Ac = mu > c

        # Obtain estimated Ac
        AcHat = muHat > c

        # -------------------------------------------------------------------
        # Muhat (interpolated) along the true Ac boundary (dict version)
        # -------------------------------------------------------------------

        #-----
        # Mu 1
        #-----

        # Muhat1 along Ac1 boundary
        muHat1_Ac1_bdry = get_bdry_values(muHat[:1,:,:], Ac1_bdry_locs)

        # Interpolate along Ac1 boundary
        muHat1_Ac1_bdry = get_bdry_vals_interpolated(muHat1_Ac1_bdry, Ac1_bdry_weights)

        #-----
        # Mu 2
        #-----

        # Muhat2 along Ac2 boundary
        muHat2_Ac2_bdry = get_bdry_values(muHat[1:,:,:], Ac2_bdry_locs)

        # Interpolate along Ac2 boundary
        muHat2_Ac2_bdry = get_bdry_vals_interpolated(muHat2_Ac2_bdry, Ac2_bdry_weights)

        #------
        # Union
        #------

        # Muhat union along Ac union boundary
        muHat_cup_Ac_cup_bdry = get_bdry_values(muHat_cup, Ac_cup_bdry_locs)

        # Interpolate along Ac union boundary
        muHat_cup_Ac_cup_bdry = get_bdry_vals_interpolated(muHat_cup_Ac_cup_bdry, Ac_cup_bdry_weights)

        #-------------
        # Intersection
        #-------------

        # Muhat intersection along Ac intersection boundary
        muHat_cap_Ac_cap_bdry = get_bdry_values(muHat_cap, Ac_cap_bdry_locs)

        # Interpolate along Ac intersection boundary
        muHat_cap_Ac_cap_bdry = get_bdry_vals_interpolated(muHat_cap_Ac_cap_bdry, Ac_cap_bdry_weights)

        # -------------------------------------------------------------------
        # Muhat (interpolated) along the true Ac boundary (concatenated
        # version)
        # -------------------------------------------------------------------

        #-----
        # Mu 1
        #-----

        # Muhat1 along Ac1 boundary
        muHat1_Ac1_bdry_concat = get_bdry_values_concat(muHat[:1,:,:], Ac1_bdry_locs)

        # Interpolate along Ac1 boundary
        muHat1_Ac1_bdry_concat = get_bdry_vals_interpolated_concat(muHat1_Ac1_bdry_concat, Ac1_bdry_weights_concat)

        #-----
        # Mu 2
        #-----

        # Muhat2 along Ac2 boundary
        muHat2_Ac2_bdry_concat = get_bdry_values_concat(muHat[1:,:,:], Ac2_bdry_locs)

        # Interpolate along Ac2 boundary
        muHat2_Ac2_bdry_concat = get_bdry_vals_interpolated_concat(muHat2_Ac2_bdry_concat, Ac2_bdry_weights_concat)

        #------
        # Union
        #------

        # Muhat union along Ac union boundary
        muHat_cup_Ac_cup_bdry_concat = get_bdry_values_concat(muHat_cup, Ac_cup_bdry_locs)

        # Interpolate along Ac union boundary
        muHat_cup_Ac_cup_bdry_concat = get_bdry_vals_interpolated_concat(muHat_cup_Ac_cup_bdry_concat, Ac_cup_bdry_weights_concat)

        #-------------
        # Intersection
        #-------------

        # Muhat intersection along Ac intersection boundary
        muHat_cap_Ac_cap_bdry_concat = get_bdry_values_concat(muHat_cap, Ac_cap_bdry_locs)

        # Interpolate along Ac intersection boundary
        muHat_cap_Ac_cap_bdry_concat = get_bdry_vals_interpolated_concat(muHat_cap_Ac_cap_bdry_concat, Ac_cap_bdry_weights_concat)

        # -------------------------------------------------------------------
        # Bootstrap 
        # -------------------------------------------------------------------
        # Initialize empty bootstrap stores for mu 1
        max_g_Ac1 = np.zeros(nBoot)
        max_g_AcHat1 = np.zeros(nBoot)

        # Initialize empty bootstrap stores for mu 2
        max_g_Ac2 = np.zeros(nBoot)
        max_g_AcHat2 = np.zeros(nBoot)

        # Initialize empty bootstrap stores for mu intersection
        max_g_Ac_cap = np.zeros(nBoot)
        max_g_AcHat_cap = np.zeros(nBoot)

        # Initialize empty bootstrap stores for mu union
        max_g_Ac_cup = np.zeros(nBoot)
        max_g_AcHat_cup = np.zeros(nBoot)

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

            # ---------------------------------------------------------------
            # Bootstrap residuals along Ac and AcHat 1
            # ---------------------------------------------------------------

            # Bootstrap residuals along Ac 1
            boot_resid1_Ac1_bdry_concat = boot_vars*resid1_Ac1_bdry_concat

            # Bootstrap residuals along AcHat 1
            boot_resid1_AcHat1_bdry_concat = boot_vars*resid1_AcHat1_bdry_concat

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of Ac1. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g_Ac1_bdry_concat = np.zeros(boot_resid1_Ac1_bdry_concat.shape[-2:])
            boot_g_Ac1_bdry_concat[...,0] = np.sum(boot_resid1_Ac1_bdry_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g_Ac1_bdry_concat[...,1] = np.sum(boot_resid1_Ac1_bdry_concat[...,1], axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along Ac1. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma1_boot_Ac1_concat = np.zeros(boot_resid1_Ac1_bdry_concat.shape[-2:])
            sigma1_boot_Ac1_concat[...,0] = np.std(boot_resid1_Ac1_bdry_concat[...,0], axis=0, ddof=1)
            sigma1_boot_Ac1_concat[...,1] = np.std(boot_resid1_Ac1_bdry_concat[...,1], axis=0, ddof=1)

            # Divide by the boostrap standard deviation on Ac
            boot_g_Ac1_bdry_concat = boot_g_Ac1_bdry_concat/sigma1_boot_Ac1_concat

            # Sum across subjects to get the bootstrapped g values along
            # the boundary of AcHat1. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g_AcHat1_bdry_concat = np.zeros(boot_resid1_AcHat1_bdry_concat.shape[-2:])
            boot_g_AcHat1_bdry_concat[...,0] = np.sum(boot_resid1_AcHat1_bdry_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g_AcHat1_bdry_concat[...,1] = np.sum(boot_resid1_AcHat1_bdry_concat[...,1], axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along AcHat1. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma1_boot_AcHat1_concat = np.zeros(boot_resid1_AcHat1_bdry_concat.shape[-2:])
            sigma1_boot_AcHat1_concat[...,0] = np.std(boot_resid1_AcHat1_bdry_concat[...,0], axis=0, ddof=1)
            sigma1_boot_AcHat1_concat[...,1] = np.std(boot_resid1_AcHat1_bdry_concat[...,1], axis=0, ddof=1)

            # Divide by the boostrap standard deviation on AcHat1
            boot_g_AcHat1_bdry_concat = boot_g_AcHat1_bdry_concat/sigma1_boot_AcHat1_concat

            # print(sigma1_boot_AcHat1_concat.shape)
            # print(sigma1_boot_AcHat1_concat)

            # print('loc', np.max(np.abs(boot_g_AcHat1_bdry_concat)))

            # Interpolation for Ac1 boundary
            boot_g_Ac1_bdry_concat = get_bdry_vals_interpolated_concat(boot_g_Ac1_bdry_concat,Ac1_bdry_weights_concat)

            # Interpolation for AcHat1 boundary
            boot_g_AcHat1_bdry_concat = get_bdry_vals_interpolated_concat(boot_g_AcHat1_bdry_concat,AcHat1_bdry_weights_concat)

            # print('mfiwaneou4an')
            # print(boot_g_AcHat1_bdry_concat.shape)
            # print(boot_g_AcHat1_bdry_concat)

            # Get maximum along Ac1 boudary
            max_g_Ac1[b] = np.max(np.abs(boot_g_Ac1_bdry_concat)) 

            # Get maximum along AcHat1 boudary
            max_g_AcHat1[b] = np.max(np.abs(boot_g_AcHat1_bdry_concat)) 

            print('!!!!!!!!!!!!!!!!!!', max_g_AcHat1[b])

            # ---------------------------------------------------------------
            # Bootstrap residuals along Ac and AcHat 2
            # ---------------------------------------------------------------

            # Bootstrap residuals along Ac 2
            boot_resid2_Ac2_bdry_concat = boot_vars*resid2_Ac2_bdry_concat

            # Bootstrap residuals along AcHat 2
            boot_resid2_AcHat2_bdry_concat = boot_vars*resid2_AcHat2_bdry_concat

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of Ac2. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g_Ac2_bdry_concat = np.zeros(boot_resid2_Ac2_bdry_concat.shape[-2:])
            boot_g_Ac2_bdry_concat[...,0] = np.sum(boot_resid2_Ac2_bdry_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g_Ac2_bdry_concat[...,1] = np.sum(boot_resid2_Ac2_bdry_concat[...,1], axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along Ac2. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma2_boot_Ac2_concat = np.zeros(boot_resid2_Ac2_bdry_concat.shape[-2:])
            sigma2_boot_Ac2_concat[...,0] = np.std(boot_resid2_Ac2_bdry_concat[...,0], axis=0, ddof=1)
            sigma2_boot_Ac2_concat[...,1] = np.std(boot_resid2_Ac2_bdry_concat[...,1], axis=0, ddof=1)

            # Divide by the boostrap standard deviation on Ac
            boot_g_Ac2_bdry_concat = boot_g_Ac2_bdry_concat/sigma2_boot_Ac2_concat

            # Sum across subjects to get the bootstrapped g values along
            # the boundary of AcHat2. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g_AcHat2_bdry_concat = np.zeros(boot_resid2_AcHat2_bdry_concat.shape[-2:])
            boot_g_AcHat2_bdry_concat[...,0] = np.sum(boot_resid2_AcHat2_bdry_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g_AcHat2_bdry_concat[...,1] = np.sum(boot_resid2_AcHat2_bdry_concat[...,1], axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along AcHat2. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma2_boot_AcHat2_concat = np.zeros(boot_resid2_AcHat2_bdry_concat.shape[-2:])
            sigma2_boot_AcHat2_concat[...,0] = np.std(boot_resid2_AcHat2_bdry_concat[...,0], axis=0, ddof=1)
            sigma2_boot_AcHat2_concat[...,1] = np.std(boot_resid2_AcHat2_bdry_concat[...,1], axis=0, ddof=1)

            # Divide by the boostrap standard deviation on AcHat2
            boot_g_AcHat2_bdry_concat = boot_g_AcHat2_bdry_concat/sigma2_boot_AcHat2_concat

            # Interpolation for Ac2 boundary
            boot_g_Ac2_bdry_concat = get_bdry_vals_interpolated_concat(boot_g_Ac2_bdry_concat,Ac2_bdry_weights_concat)

            # Interpolation for AcHat2 boundary
            boot_g_AcHat2_bdry_concat = get_bdry_vals_interpolated_concat(boot_g_AcHat2_bdry_concat,AcHat2_bdry_weights_concat)

            # Get maximum along Ac2 boudary
            max_g_Ac2[b] = np.max(np.abs(boot_g_Ac2_bdry_concat)) 

            # Get maximum along AcHat2 boudary
            max_g_AcHat2[b] = np.max(np.abs(boot_g_AcHat2_bdry_concat)) 

            # ---------------------------------------------------------------
            # Bootstrap residuals along Ac intersection 
            # ---------------------------------------------------------------
            # Bootstrap residuals along Ac intersection.
            boot_resid_Ac_cap_bdry_concat = boot_vars.reshape((*boot_vars.shape),1)*resid_Ac_cap_bdry_concat

            # Bootstrap residuals along AcHat intersection.
            boot_resid_AcHat_cap_bdry_concat = boot_vars.reshape((*boot_vars.shape),1)*resid_AcHat_cap_bdry_concat

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of Ac intersection. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g_Ac_cap_bdry_concat = np.zeros(boot_resid_Ac_cap_bdry_concat.shape[-3:])
            boot_g_Ac_cap_bdry_concat[...,0] = np.sum(boot_resid_Ac_cap_bdry_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g_Ac_cap_bdry_concat[...,1] = np.sum(boot_resid_Ac_cap_bdry_concat[...,1], axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along Ac_cap. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma_boot_Ac_cap_concat = np.zeros(boot_resid_Ac_cap_bdry_concat.shape[-3:])
            sigma_boot_Ac_cap_concat[...,0] = np.std(boot_resid_Ac_cap_bdry_concat[...,0], axis=0, ddof=1)
            sigma_boot_Ac_cap_concat[...,1] = np.std(boot_resid_Ac_cap_bdry_concat[...,1], axis=0, ddof=1)

            # Divide by the boostrap standard deviation on Ac
            boot_g_Ac_cap_bdry_concat = boot_g_Ac_cap_bdry_concat/sigma_boot_Ac_cap_concat

            # Sum across subjects to get the bootstrapped g values along
            # the boundary of AcHat intersection. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g_AcHat_cap_bdry_concat = np.zeros(boot_resid_AcHat_cap_bdry_concat.shape[-3:])
            boot_g_AcHat_cap_bdry_concat[...,0] = np.sum(boot_resid_AcHat_cap_bdry_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g_AcHat_cap_bdry_concat[...,1] = np.sum(boot_resid_AcHat_cap_bdry_concat[...,1], axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along AcHat intersection. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma_boot_AcHat_cap_concat = np.zeros(boot_resid_AcHat_cap_bdry_concat.shape[-3:])
            sigma_boot_AcHat_cap_concat[...,0] = np.std(boot_resid_AcHat_cap_bdry_concat[...,0], axis=0, ddof=1)
            sigma_boot_AcHat_cap_concat[...,1] = np.std(boot_resid_AcHat_cap_bdry_concat[...,1], axis=0, ddof=1)

            # Divide by the boostrap standard deviation on AcHat intersection.
            boot_g_AcHat_cap_bdry_concat = boot_g_AcHat_cap_bdry_concat/sigma_boot_AcHat_cap_concat

            # Interpolation for Ac intersection boundary
            boot_g_Ac_cap_bdry_concat = get_bdry_vals_interpolated_concat(boot_g_Ac_cap_bdry_concat,Ac_cap_bdry_weights_concat)

            # Interpolation for AcHat intersection boundary
            boot_g_AcHat_cap_bdry_concat = get_bdry_vals_interpolated_concat(boot_g_AcHat_cap_bdry_concat,AcHat_cap_bdry_weights_concat)

            # Get maximum along Ac intersection boudary
            max_g_Ac_cap[b] = np.max(np.abs(boot_g_Ac_cap_bdry_concat)) 

            # Get maximum along AcHat intersection boudary
            max_g_AcHat_cap[b] = np.max(np.abs(boot_g_AcHat_cap_bdry_concat)) 

            # ---------------------------------------------------------------
            # Bootstrap residuals along Ac union 
            # ---------------------------------------------------------------

            # Bootstrap residuals along Ac union.
            boot_resid_Ac_cup_bdry_concat = boot_vars.reshape((*boot_vars.shape),1)*resid_Ac_cup_bdry_concat

            # Bootstrap residuals along AcHat union.
            boot_resid_AcHat_cup_bdry_concat = boot_vars.reshape((*boot_vars.shape),1)*resid_AcHat_cup_bdry_concat

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of Ac union. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g_Ac_cup_bdry_concat = np.zeros(boot_resid_Ac_cup_bdry_concat.shape[-3:])
            boot_g_Ac_cup_bdry_concat[...,0] = np.sum(boot_resid_Ac_cup_bdry_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g_Ac_cup_bdry_concat[...,1] = np.sum(boot_resid_Ac_cup_bdry_concat[...,1], axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along Ac_cup. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma_boot_Ac_cup_concat = np.zeros(boot_resid_Ac_cup_bdry_concat.shape[-3:])
            sigma_boot_Ac_cup_concat[...,0] = np.std(boot_resid_Ac_cup_bdry_concat[...,0], axis=0, ddof=1)
            sigma_boot_Ac_cup_concat[...,1] = np.std(boot_resid_Ac_cup_bdry_concat[...,1], axis=0, ddof=1)

            # Divide by the boostrap standard deviation on Ac
            boot_g_Ac_cup_bdry_concat = boot_g_Ac_cup_bdry_concat/sigma_boot_Ac_cup_concat

            # Sum across subjects to get the bootstrapped g values along
            # the boundary of AcHat union. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g_AcHat_cup_bdry_concat = np.zeros(boot_resid_AcHat_cup_bdry_concat.shape[-3:])
            boot_g_AcHat_cup_bdry_concat[...,0] = np.sum(boot_resid_AcHat_cup_bdry_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g_AcHat_cup_bdry_concat[...,1] = np.sum(boot_resid_AcHat_cup_bdry_concat[...,1], axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along AcHat union. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma_boot_AcHat_cup_concat = np.zeros(boot_resid_AcHat_cup_bdry_concat.shape[-3:])
            sigma_boot_AcHat_cup_concat[...,0] = np.std(boot_resid_AcHat_cup_bdry_concat[...,0], axis=0, ddof=1)
            sigma_boot_AcHat_cup_concat[...,1] = np.std(boot_resid_AcHat_cup_bdry_concat[...,1], axis=0, ddof=1)

            # Divide by the boostrap standard deviation on AcHat union.
            boot_g_AcHat_cup_bdry_concat = boot_g_AcHat_cup_bdry_concat/sigma_boot_AcHat_cup_concat

            # Interpolation for Ac union boundary
            boot_g_Ac_cup_bdry_concat = get_bdry_vals_interpolated_concat(boot_g_Ac_cup_bdry_concat,Ac_cup_bdry_weights_concat)

            # Interpolation for AcHat union boundary
            boot_g_AcHat_cup_bdry_concat = get_bdry_vals_interpolated_concat(boot_g_AcHat_cup_bdry_concat,AcHat_cup_bdry_weights_concat)

            # Get maximum along Ac union boudary
            max_g_Ac_cup[b] = np.max(np.abs(boot_g_Ac_cup_bdry_concat))  

            # Get maximum along AcHat union boudary
            max_g_AcHat_cup[b] = np.max(np.abs(boot_g_AcHat_cup_bdry_concat)) 


        # Get the elementwise maximum of the a's from the mu1 and mu2
        # boundaries for each realization.
        max_g_Ac = np.maximum(max_g_Ac1,max_g_Ac2) 

        # Get maximum along AcHat2 boudary
        max_g_AcHat = np.maximum(max_g_AcHat1,max_g_AcHat2) 

        print('here2')
        print(max_g_AcHat)
        print(max_g_AcHat1,max_g_AcHat2)

        t2 = time.time()
        print('Bootstrap time: ', t2-t1)

        # -------------------------------------------------------------------
        # Get the statistic field which defined FcHat^{+/-}=AcHat1^{+/-}
        # intersect AcHat2^{+/-}
        # -------------------------------------------------------------------
        # Get the statistic field which defined Achat^{+/-} for intersection
        stat = ((muHat-c)/(sigma*tau)).reshape(1,(*muHat.shape))

        print('stat shape', stat.shape)
        stat = np.min(stat, axis=1)

        print('sig shape', (sigma).shape)
        print('tau shape', (tau).shape)
        print('sigtau shape', (sigma*tau).shape)

        print('stat shape', stat.shape)

        # -------------------------------------------------------------------
        # Get intersection field, Fc
        # -------------------------------------------------------------------
        Fc = np.amax(Ac,axis=0)

        # -------------------------------------------------------------------
        # Obtaining a from percentiles of the max distribution for doing the
        # two boundaries seperately
        # -------------------------------------------------------------------

        # Get the a estimates along the true Ac boundary and estimated AcHat boundary
        a_trueBdry = np.percentile(max_g_Ac, 100*p).reshape(nPvals,1,1,1)
        a_estBdry = np.percentile(max_g_AcHat, 100*p).reshape(nPvals,1,1,1) # [pvals, 1, [1 for _ in dim>1]]

        # Reformat them to an array form useful for boolean operation
        a_trueBdry = np.concatenate((-a_trueBdry,a_trueBdry),axis=1)
        a_estBdry = np.concatenate((-a_estBdry,a_estBdry),axis=1)

        # Obtain AcHat^+ and AcHat^- based on a from the true boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        FcHat_pm_trueBdry = stat >= a_trueBdry

        # Obtain AcHat^+ and AcHat^- based on a from the estimated boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        FcHat_pm_estBdry = stat >= a_estBdry

        # -------------------------------------------------------------------
        # Obtaining a from percentiles of the max distribution for doing the
        # union boundary
        # -------------------------------------------------------------------

        # Get the a estimates along the true Ac boundary and estimated AcHat boundary
        a_trueBdry_cup = np.percentile(max_g_Ac_cup, 100*p).reshape(nPvals,1,1,1)
        a_estBdry_cup = np.percentile(max_g_AcHat_cup, 100*p).reshape(nPvals,1,1,1) # [pvals, 1, [1 for _ in dim>1]]

        # Reformat them to an array form useful for boolean operation
        a_trueBdry_cup = np.concatenate((-a_trueBdry_cup,a_trueBdry_cup),axis=1)
        a_estBdry_cup = np.concatenate((-a_estBdry_cup,a_estBdry_cup),axis=1)

        print('wefuqqubf')
        print(stat.shape)
        print(a_trueBdry_cup.shape)

        # Obtain AcHat^+ and AcHat^- based on a from the true boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        FcHat_pm_trueBdry_cup = stat >= a_trueBdry_cup

        print(FcHat_pm_trueBdry_cup.shape)

        # Obtain AcHat^+ and AcHat^- based on a from the estimated boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        FcHat_pm_estBdry_cup = stat >= a_estBdry_cup

        # -------------------------------------------------------------------
        # Obtaining a from percentiles of the max distribution for doing the
        # intersection boundary
        # -------------------------------------------------------------------

        # Get the a estimates along the true boundary and estimated boundary
        a_trueBdry_cap = np.percentile(max_g_Ac_cap, 100*p).reshape(nPvals,1,1,1)
        a_estBdry_cap = np.percentile(max_g_AcHat_cap, 100*p).reshape(nPvals,1,1,1) # [pvals, 1, [1 for _ in dim>1]]

        # Reformat them to an array form useful for boolean operation
        a_trueBdry_cap = np.concatenate((-a_trueBdry_cap,a_trueBdry_cap),axis=1)
        a_estBdry_cap = np.concatenate((-a_estBdry_cap,a_estBdry_cap),axis=1)

        # Obtain AcHat^+ and AcHat^- based on a from the true boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        FcHat_pm_trueBdry_cap = stat >= a_trueBdry_cap

        # Obtain AcHat^+ and AcHat^- based on a from the estimated boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        FcHat_pm_estBdry_cap = stat >= a_estBdry_cap

        # -------------------------------------------------------------------
        # Some set logic to work out violations for combination of boundaries
        # estimated a
        # -------------------------------------------------------------------

        # Obtain AcHat^+\AcHat based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        FcHatp_sub_Fc_trueBdry = FcHat_pm_trueBdry[:,1,...] & ~Fc

        # Obtain AcHat^+\AcHat based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        FcHatp_sub_Fc_estBdry = FcHat_pm_estBdry[:,1,...] & ~Fc

        # Obtain Ac\AcHat^- based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        Fc_sub_FcHatm_trueBdry = Fc & ~FcHat_pm_trueBdry[:,0,...]

        # Obtain Ac\AcHat^- based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        Fc_sub_FcHatm_estBdry = Fc & ~FcHat_pm_estBdry[:,0,...]

        # -------------------------------------------------------------------
        # Some set logic to work out violations for union boundary estimated 
        # a
        # -------------------------------------------------------------------

        # Obtain AcHat^+\AcHat based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        FcHatp_sub_Fc_trueBdry_cup = FcHat_pm_trueBdry_cup[:,1,...] & ~Fc

        # Obtain AcHat^+\AcHat based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        FcHatp_sub_Fc_estBdry_cup = FcHat_pm_estBdry_cup[:,1,...] & ~Fc

        # Obtain Ac\AcHat^- based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        Fc_sub_FcHatm_trueBdry_cup = Fc & ~FcHat_pm_trueBdry_cup[:,0,...]

        # Obtain Ac\AcHat^- based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        Fc_sub_FcHatm_estBdry_cup = Fc & ~FcHat_pm_estBdry_cup[:,0,...]

        # -------------------------------------------------------------------
        # Some set logic to work out violations for intersection boundary
        # estimated a
        # -------------------------------------------------------------------

        # Obtain AcHat^+\AcHat based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        FcHatp_sub_Fc_trueBdry_cap = FcHat_pm_trueBdry_cap[:,1,...] & ~Fc

        # Obtain AcHat^+\AcHat based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        FcHatp_sub_Fc_estBdry_cap = FcHat_pm_estBdry_cap[:,1,...] & ~Fc

        # Obtain Ac\AcHat^- based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        Fc_sub_FcHatm_trueBdry_cap = Fc & ~FcHat_pm_trueBdry_cap[:,0,...]

        # Obtain Ac\AcHat^- based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        Fc_sub_FcHatm_estBdry_cap = Fc & ~FcHat_pm_estBdry_cap[:,0,...]

        # -------------------------------------------------------------------
        # Check whether there were any boundary violations using voxelwise
        # set logic (checking if voxels existed in one set but not another, 
        # etc) for the seperate boundary version.
        # -------------------------------------------------------------------

        # Record if we saw a violation in the true boundary based sets
        trueBdry_success[r,:] = 1-(np.any(FcHatp_sub_Fc_trueBdry,axis=(1,2)) | np.any(Fc_sub_FcHatm_trueBdry,axis=(1,2))) # MARKER: AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success[r,:] = 1-(np.any(FcHatp_sub_Fc_estBdry,axis=(1,2)) | np.any(Fc_sub_FcHatm_estBdry,axis=(1,2)))# MARKER: AXES WONT WORK FOR 3D ATM

        # -------------------------------------------------------------------
        # Check whether there were any boundary violations using voxelwise
        # set logic (checking if voxels existed in one set but not another, 
        # etc) for the union boundary version.
        # -------------------------------------------------------------------

        # Record if we saw a violation in the true boundary based sets
        trueBdry_success_cup[r,:] = 1-(np.any(FcHatp_sub_Fc_trueBdry_cup,axis=(1,2)) | np.any(Fc_sub_FcHatm_trueBdry_cup,axis=(1,2))) # MARKER: AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success_cup[r,:] = 1-(np.any(FcHatp_sub_Fc_estBdry_cup,axis=(1,2)) | np.any(Fc_sub_FcHatm_estBdry_cup,axis=(1,2)))# MARKER: AXES WONT WORK FOR 3D ATM

        # -------------------------------------------------------------------
        # Check whether there were any boundary violations using voxelwise
        # set logic (checking if voxels existed in one set but not another, 
        # etc) for the intersection boundary version.
        # -------------------------------------------------------------------

        # Record if we saw a violation in the true boundary based sets
        trueBdry_success_cap[r,:] = 1-(np.any(FcHatp_sub_Fc_trueBdry_cap,axis=(1,2)) | np.any(Fc_sub_FcHatm_trueBdry_cap,axis=(1,2))) # MARKER: AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success_cap[r,:] = 1-(np.any(FcHatp_sub_Fc_estBdry_cap,axis=(1,2)) | np.any(Fc_sub_FcHatm_estBdry_cap,axis=(1,2)))# MARKER: AXES WONT WORK FOR 3D ATM


        # -------------------------------------------------------------------
        # Get stat along the Fc boundary
        # -------------------------------------------------------------------
        # Get the values along the outer and inner boundaries
        stat_FcBdry = get_bdry_values(stat, Ac_cap_bdry_locs)

        # Interpolate to get the values along the true boundary
        stat_FcBdry = get_bdry_vals_interpolated(stat_FcBdry, Fc_bdry_weights)

        print('F: ', stat_FcBdry)
        print('est bdry 1: ', a_estBdry[:,0,:,0])
        print('est bdry 2: ', a_estBdry[:,1,:,0])
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

        # -------------------------------------------------------------------
        # Check whether there were any boundary violations using interpolated
        # boundary values (checking if voxels had values corresponding to no
        # violations, etc) for union boundary
        # -------------------------------------------------------------------

        # Perform lower check on stat map using thresholds based on the
        # estimated boundary
        bdry_lowerCheck_estBdry_cup = stat_FcBdry >= a_estBdry_cup[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # estimated boundary
        bdry_upperCheck_estBdry_cup = stat_FcBdry <= a_estBdry_cup[:,1,:,0]

        # Perform lower check on stat map using thresholds based on the
        # true boundary
        bdry_lowerCheck_trueBdry_cup = stat_FcBdry >= a_trueBdry_cup[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # true boundary
        bdry_upperCheck_trueBdry_cup = stat_FcBdry <= a_trueBdry_cup[:,1,:,0]

        # -------------------------------------------------------------------
        # Check whether there were any boundary violations using interpolated
        # boundary values (checking if voxels had values corresponding to no
        # violations, etc) for intersection boundary
        # -------------------------------------------------------------------

        # Perform lower check on stat map using thresholds based on the
        # estimated boundary
        bdry_lowerCheck_estBdry_cap = stat_FcBdry >= a_estBdry_cap[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # estimated boundary
        bdry_upperCheck_estBdry_cap = stat_FcBdry <= a_estBdry_cap[:,1,:,0]

        # Perform lower check on stat map using thresholds based on the
        # true boundary
        bdry_lowerCheck_trueBdry_cap = stat_FcBdry >= a_trueBdry_cap[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # true boundary
        bdry_upperCheck_trueBdry_cap = stat_FcBdry <= a_trueBdry_cap[:,1,:,0]

        # -------------------------------------------------------------------
        # Work out whether simulation observed successful sets.
        # -------------------------------------------------------------------
        # Record if we saw a violation in the true boundary based sets
        trueBdry_success_intrp[r,:] = np.all(bdry_lowerCheck_trueBdry,axis=(1)) & np.all(bdry_upperCheck_trueBdry,axis=(1))# MARKER: AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success_intrp[r,:] = np.all(bdry_lowerCheck_estBdry,axis=(1)) & np.all(bdry_upperCheck_estBdry,axis=(1)) # MARKER: AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the true boundary based sets
        trueBdry_success_intrp_cup[r,:] = np.all(bdry_lowerCheck_trueBdry_cup,axis=(1)) & np.all(bdry_upperCheck_trueBdry_cup,axis=(1))# MARKER: AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success_intrp_cup[r,:] = np.all(bdry_lowerCheck_estBdry_cup,axis=(1)) & np.all(bdry_upperCheck_estBdry_cup,axis=(1)) # MARKER: AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the true boundary based sets
        trueBdry_success_intrp_cap[r,:] = np.all(bdry_lowerCheck_trueBdry_cap,axis=(1)) & np.all(bdry_upperCheck_trueBdry_cap,axis=(1))# MARKER: AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success_intrp_cap[r,:] = np.all(bdry_lowerCheck_estBdry_cap,axis=(1)) & np.all(bdry_upperCheck_estBdry_cap,axis=(1)) # MARKER: AXES WONT WORK FOR 3D ATM


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# HERE

    # For the interpolated boundary success checks, we still need to do the 
    # voxelwise checks as well. This will take care of that.
    trueBdry_success_intrp = trueBdry_success_intrp*trueBdry_success
    estBdry_success_intrp = estBdry_success_intrp*estBdry_success

    # For the interpolated boundary success checks, we still need to do the 
    # voxelwise checks as well. This will take care of that.
    trueBdry_success_intrp_cup = trueBdry_success_intrp_cup*trueBdry_success_cup
    estBdry_success_intrp_cup = estBdry_success_intrp_cup*estBdry_success_cup

    # For the interpolated boundary success checks, we still need to do the 
    # voxelwise checks as well. This will take care of that.
    trueBdry_success_intrp_cap = trueBdry_success_intrp_cap*trueBdry_success_cap
    estBdry_success_intrp_cap = estBdry_success_intrp_cap*estBdry_success_cap

    # Coverage probabilities
    coverage_trueBdry = np.mean(trueBdry_success,axis=0)
    coverage_estBdry = np.mean(estBdry_success,axis=0)

    # Coverage probabilities
    coverage_trueBdry_cup = np.mean(trueBdry_success_cup,axis=0)
    coverage_estBdry_cup = np.mean(estBdry_success_cup,axis=0)

    # Coverage probabilities
    coverage_trueBdry_cap = np.mean(trueBdry_success_cap,axis=0)
    coverage_estBdry_cap = np.mean(estBdry_success_cap,axis=0)

    # Coverage probabilities
    coverage_trueBdry_intrp = np.mean(trueBdry_success_intrp,axis=0)
    coverage_estBdry_intrp = np.mean(estBdry_success_intrp,axis=0)

    # Coverage probabilities
    coverage_trueBdry_intrp_cup = np.mean(trueBdry_success_intrp_cup,axis=0)
    coverage_estBdry_intrp_cup = np.mean(estBdry_success_intrp_cup,axis=0)

    # Coverage probabilities
    coverage_trueBdry_intrp_cap = np.mean(trueBdry_success_intrp_cap,axis=0)
    coverage_estBdry_intrp_cap = np.mean(estBdry_success_intrp_cap,axis=0)

    print('Coverage: ', coverage_trueBdry_intrp)
    print('Coverage_cup: ', coverage_trueBdry_intrp_cup)
    print('Coverage_cap: ', coverage_trueBdry_intrp_cap)

    # Save the violations to a file
    append_to_file('trueSuccess'+str(nSub)+'.csv', trueBdry_success) 
    append_to_file('estSuccess'+str(nSub)+'.csv', estBdry_success)
    append_to_file('trueSuccess'+str(nSub)+'_intrp.csv', trueBdry_success_intrp) 
    append_to_file('estSuccess'+str(nSub)+'_intrp.csv', estBdry_success_intrp)

    t2overall = time.time()

    print('overall time: ', t2overall-t1overall)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Run example
SpatialSims2mu('/home/tommaullin/Documents/ConfSets/',100, {'type': 'circle2D_2mu', 'center': np.array([[-20,0],[20,0]]), 'fwhm': np.array([3,3]), 'r': np.array([25,25]), 'mag': np.array([3,3])}, 5, 2, np.linspace(0,1,21))
#SpatialSims('/home/tommaullin/Documents/ConfSets/',100, {'type': 'circle2D', 'center': np.array([0,0]), 'fwhm': np.array([5,5]), 'r': 30, 'mag': 3}, 1, 2, np.linspace(0,1,21))

#SpatialSims('/home/tommaullin/Documents/ConfSets/',100, {'type': 'circle2D', 'center': np.array([0,0]), 'fwhm': np.array([3,3]), 'r': 30, 'mag': 3}, 200, 2, np.array([0.8,0.9,0.95]))