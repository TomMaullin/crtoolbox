import os
import sys
import time
import numpy as np
from lib.boundary import *
from lib.fileio import *
from bootstrap import *
from lib.set_theory import powerset
import yaml
import matplotlib.pyplot as plt


"""
This function generates confidence regions for conjunction inference. 

Inputs:
    betahat_files: A list of strings representing the paths to the betahat files.
    sigmahat_file: A string representing the path to the sigmahat file.
    resid_files: A list of strings representing the paths to the residual files.
    out_dir: A string representing the path to the output directory.
    c: A float representing the threshold for the conjunction inference.
    p: An array of floats representing the desired coverage of the confidence
       regions.
    m: An integer representing the number of fields we're estimating the conjunction
       of.
    mask: A numpy array of shape (x,y,z) representing the mask for the data. This
          is not currently implemented.
    n_boot: An integer representing the number of bootstrap samples to use.
    tau: A string representing the value of tau to use as a function of n_sub; it
         defaults to 1/sqrt(n_sub).
    output: A boolean representing whether to output the confidence regions as nifti

Outputs:
    FcHat_plus_files: A list of strings representing the paths to the upper confidence
                      regions.
    FcHat_minus_files: A list of strings representing the paths to the lower confidence 
                       regions.
    FcHat_files: A string representing the path to the estimated conjunction region.
    a_estBdry: A numpy array of shape (n_boot,) representing the bootstrap quantiles.
""" 
def generate_CRs(mean_fname, sig_fname, res_fnames, out_dir, c, p, m=1, mask=None, n_boot=5000, tau='1/np.sqrt(n_sub)', output=True):

    # If p is not an array, make it one
    if not isinstance(p, np.ndarray):    
        p = np.array([p])

    # Get number of subjects
    n_sub = len(res_fnames)

    # Read in mean and sigma
    muHats = cycle_axes(read_images(mean_fname))
    sigmas = cycle_axes(read_images(sig_fname))

    # Get image dimensions
    image_dim = muHats.shape[1:]
    D = len(image_dim)

    # Evaulate tau
    tau = eval(tau)

    # If there is no mask, make one from the zero set in sigma
    if mask is None:

        # Loop through m
        for i in np.arange(m):

            # Get sigma for this field
            sigma = sigmas[i,...]

            # Get mask
            mask_current = sigma != 0

            # If this is the first iteration, initialize mask
            if i == 0:

                mask = mask_current

            # Otherwise, take the intersection of the current mask and the previous mask
            else:

                mask = np.logical_and(mask, mask_current)
    

    # Boundary locations and values
    # Make a structure to hold the estimated boundary weights in array form 
    est_bdry_weights_concat = {}

    # Make a structure to hold the true and estimated boundary locations
    est_bdry_locs = {}

    for i in np.arange(m):

        # Get muhat for this sample
        muHat = muHats[i,...]

        # Boundary locations for AcHati
        # Get boolean maps for the boundary of AcHat
        AcHat_bdry_map = get_bdry_maps(muHat, c, mask)

        # Get coordinates for the boundary of AcHat
        AcHat_bdry_locs = get_bdry_locs(AcHat_bdry_map)

        # Save boundary locations
        est_bdry_locs['AcHat'+str(i+1)] = AcHat_bdry_locs

        # Delete map as we no longer need it
        del AcHat_bdry_map

        # Interpolation weights for AcHati boundary (Array version)
        # Obtain the values along the boundary for AcHati
        AcHat_bdry_vals_concat = get_bdry_values_concat(muHat, AcHat_bdry_locs)

        # Obtain the weights along the boundary for AcHati
        AcHat_bdry_weights_concat = get_bdry_weights_concat(AcHat_bdry_vals_concat, c)

        # Save boundary weights
        est_bdry_weights_concat['AcHat'+str(i+1)] = AcHat_bdry_weights_concat

        # Delete values as we no longer need them
        del AcHat_bdry_vals_concat


    # Get minimum fields
    # This is named cap as the excursion set of the minimum field is
    # the intersection of all fields (\cap in latex)
    cap_muHat = np.amin(muHats,axis=0)

    # Boundary locations for FcHat
    # Get boolean map for the boundary of FcHat
    FcHat_bdry_map = get_bdry_maps(cap_muHat, c, mask)

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

        # Loop through subjects
        for j in np.arange(n_sub):
        
            # Obtain residuals
            resid = cycle_axes(read_images(res_fnames[j]))

            # Residuals along FcHat boundary for current subject
            current_resid_dFcHat_concat = get_bdry_values_concat(resid, FcHat_bdry_locs)

            # Concatenate residuals
            if j == 0:
                resid_dFcHat_concat = current_resid_dFcHat_concat
            else:
                resid_dFcHat_concat = np.concatenate((resid_dFcHat_concat, current_resid_dFcHat_concat), axis=0)

        # Save residuals
        resids_dFcHat['field'+str(i+1)] = resid_dFcHat_concat

        # Obtain MuHat along FcHat
        muHat_dFcHat_concat = get_bdry_values_concat(muHats[i,...], FcHat_bdry_locs)

        # Save mu
        muHat_dFcHat['field'+str(i+1)] = muHat_dFcHat_concat

    # Delete residual as it is longer needed
    del resid

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

    # Obtain FcHat
    FcHat = cap_muHat > c

    # Perform Bootstrap 
    a_estBdry = bootstrap_resids(resids_dFcHat_partitioned, weights_dFcHat, m, n_boot, p, n_sub)

    # Reshape a_estBdry to be the same dimensions as before, followed by 1 
    # for each dimension of the field
    a_estBdry = a_estBdry.reshape(a_estBdry.shape+tuple(np.ones(D,dtype=int)))

    # Create empty g array
    g = np.zeros(muHats.shape)

    # Reform mask to be the same shape as muHats
    mask = np.broadcast_to(mask, muHats.shape)

    # Get the statistic field which defined Achat^{+/-,i}
    g[mask] = ((muHats[mask]-c)/(sigmas[mask]*tau))

    # Take minimum over i
    stat = np.amin(g,axis=0)

    # If we are outputting CRs
    if output:

        # Obtain FcHat^+ and FcHat^- based on a from the estimated boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        FcHat_pm_estBdry = stat >= a_estBdry

        # Apply mask
        FcHat_pm_estBdry = FcHat_pm_estBdry*mask

        # Empty stores for filenames
        FcHat_plus_files = []
        FcHat_minus_files = []
        FcHat_files = []

        # Loop through alpha values saving confidence regions
        for i in np.arange(len(p)):

            # If D is 3, output nifti image
            if D == 3:

                # Save upper confidence region
                addBlockToNifti(os.path.join(out_dir,"Upper_CR_"+str(p[i])+".nii"),
                                FcHat_pm_estBdry[i,0,...], np.arange(np.prod(image_dim)),
                                volInd=0,dim=image_dim)
                
                # Save filename
                FcHat_plus_files.append(os.path.join(out_dir,"Upper_CR_"+str(p[i])+".nii"))
                
                # Save lower confidence region
                addBlockToNifti(os.path.join(out_dir,"Lower_CR_"+str(p[i])+".nii"),
                                FcHat_pm_estBdry[i,1,...], np.arange(np.prod(image_dim)),
                                volInd=0,dim=image_dim)
                
                # Save filename
                FcHat_minus_files.append(os.path.join(out_dir,"Lower_CR_"+str(p[i])+".nii"))

            # Otherwise output as a numpy array
            else:

                # Save upper confidence region
                np.save(os.path.join(out_dir,"Upper_CR_"+str(p[i])+".npy"),
                        FcHat_pm_estBdry[i,0,...])
                
                # Save filename
                FcHat_plus_files.append(os.path.join(out_dir,"Upper_CR_"+str(p[i])+".npy"))

                # Save lower confidence region
                np.save(os.path.join(out_dir,"Lower_CR_"+str(p[i])+".npy"),
                        FcHat_pm_estBdry[i,1,...])
                
                # Save filename
                FcHat_minus_files.append(os.path.join(out_dir,"Lower_CR_"+str(p[i])+".npy"))
            
        # If D is 3, output nifti image
        if D == 3:
                
            # Save estimated conjunction region
            addBlockToNifti(os.path.join(out_dir,"Estimated_Ac.nii"),
                            FcHat, np.arange(np.prod(image_dim)),
                            volInd=0,dim=image_dim)
            
            # Save filename
            FcHat_files.append(os.path.join(out_dir,"Estimated_Ac.nii"))
        
        # Otherwise output as a numpy array
        else:

            # Save estimated conjunction region
            np.save(os.path.join(out_dir,"Estimated_Ac.npy"),
                    FcHat)
            
            # Save filename
            FcHat_files.append(os.path.join(out_dir,"Estimated_Ac.npy"))

        # Check if filenames are of length 1
        if len(FcHat_plus_files) == 1:
            FcHat_plus_files = FcHat_plus_files[0]
        if len(FcHat_minus_files) == 1:
            FcHat_minus_files = FcHat_minus_files[0]
        if len(FcHat_files) == 1:
            FcHat_files = FcHat_files[0]

    # Otherwise, set everything to None
    else:   
        FcHat_plus_files = None
        FcHat_minus_files = None
        FcHat_files = None

    # Reshape a_estBdry to be the same dimensions as before, followed by 1
    a_estBdry = a_estBdry[:,1,...].reshape(np.prod(a_estBdry[:,1,...].shape))

    # Return result
    return(FcHat_plus_files, FcHat_minus_files, FcHat_files, a_estBdry)
        

# Small helper function to cycle axes
def cycle_axes(arr):
    """
    Cycles the axes of an array so that the last axis becomes the first axis.
    
    Parameters
    ----------
    arr : array_like
        Array to have axes cycled.
    
    Returns
    -------
    arr : array_like
        Array with axes cycled.
    """

    # Get number of dimensions
    ndims = len(arr.shape)

    # Return array with axes cycled
    return arr.transpose(np.roll(np.arange(ndims), 1))