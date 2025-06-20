import os
import sys
import time
import numpy as np
import warnings
from crtoolbox.lib.boundary import *
from crtoolbox.lib.fileio import *
from crtoolbox.bootstrap import *
from crtoolbox.lib.set_theory import powerset
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
    tau: Image or scalar of tau_n values
    X: a numpy array of shape (n,p) representing the design matrix in a regression 
       context. If this is None then a vector of ones is used, and a simple mean
       model is assumed.
    L: a numpy array of shape (p,1) representing the contrast matrix in a regression
       model. If this is None then a (1,1) matrix with a single entry of 1 is used,
       and a simple mean model is assumed.
    output: A boolean representing whether to output the confidence regions as nifti
            files.
    mode: A string representing the type of CRs to generate. This can be 'sss' for the
          method of Sommerfeld, Sain and Schwartzman (2018), or 'cohens' for the cohens
          d method of Bowring et al. (2019). 
    method: An integer representing the method to use for the cohens d method. This can
            be 1, 2 or 3. The three algoirthms are those listed in the 
            following manuscript:
                https://doi.org/10.1016/j.neuroimage.2020.117477
            (c.f. Section 2.6)
            Method 2 is used by default following the recommendations of the above
            reference. This argument is not needed for 'sss' mode.
            Note: Method 3 is not currently implemented.
    post_fix: A string representing the post-fix to add to the output files. 
              Defaults to None, which means no post-fix is added.

Outputs:
    FcHat_plus_files: A list of strings representing the paths to the upper confidence
                      regions.
    FcHat_minus_files: A list of strings representing the paths to the lower confidence 
                       regions.
    FcHat_files: A string representing the path to the estimated conjunction region.
    a_estBdry: A numpy array of shape (n_boot,) representing the bootstrap quantiles.
""" 
def generate_CRs(mean_fname, sig_fname, res_fnames, out_dir=None, c=None, 
                 p=0.95, mask=None, n_boot=5000, tau=None, X=None, L=None,
                 output=True, mode='sss', method=2, post_fix=None):

    # Check if the mode is valid
    if mode not in ['sss', 'cohens']:
        raise ValueError('Invalid mode. Must be one of sss or cohens.')
    
    # Check if the post_fix starts with an underscore
    if post_fix is not None and not post_fix.startswith('_'):
        # If not, add an underscore to the start of the post_fix
        post_fix = '_' + post_fix
    elif post_fix is None:
        # If post_fix is None, set it to an empty string
        post_fix = ''
    
    # Check if mean_fname is a list
    if isinstance(mean_fname, list):
        # Set m equal to the length of mean_fname
        m = len(mean_fname)
    else:   
        # Set m equal to 1
        m = 1

    # If p is not an array, make it one
    if not isinstance(p, np.ndarray):    
        p = np.array([p])

    # res_fnames must be a list of lists of filenames
    if not isinstance(res_fnames, list):
        raise TypeError('residuals must be a list of filenames or a list of lists of filenames.')
    elif isinstance(res_fnames[0], str):
        res_fnames = [res_fnames for _ in range(m)]


    # Create empty subject number dict
    n_sub = {}

    # Get number of subjects
    for i in np.arange(m):
        # Get number of subjects for this condition
        n_sub[str(i)] = len(res_fnames[i])

    # Read in mean and sigma
    muHats = cycle_axes(read_images(mean_fname))

    # If the mode is cohens, adjust c
    if mode == 'cohens':

        # Check is m is greater than 1
        if m > 1:
            raise ValueError('Conjunctions for Cohens d are not currently supported.')

        # Adjust c
        c = 1/(1-3/(4*n_sub[str(0)] - 5)) * c

        # If the method is 1, set sigma to the correct value
        if method == 1:

            # Set sigma to the square root of the 1 plus d squared over 2
            sigma = np.sqrt(1 + muHats**2/2)

        # If the method is 2, check the sigma file was provided
        elif method == 2:

            # Check sigma file was provided
            if sig_fname is None:

                raise ValueError('Sigma file must be provided for Cohens d method 2.')
            
            else:
            
                sigmas = cycle_axes(read_images(sig_fname))
            
        # If the method is 3, raise an error as it is not currently implemented
        elif method == 3:

            raise NotImplementedError('Method 3 for Cohens d is not implemented.')
            
    # If the mode is sss, check sigma file was provided
    elif mode == 'sss':

        # Check sigma file was provided
        if sig_fname is None:

            raise ValueError('Sigma file must be provided for SSS method.')
    
        else:
        
            sigmas = cycle_axes(read_images(sig_fname))

    # Get image dimensions
    image_dim = muHats.shape[1:]
    D = len(image_dim)


    # Check if the mask is a filename or a numpy array.
    if isinstance(mask, str) and os.path.isfile(mask):

        # If mask is a filename, read in the file
        mask = (cycle_axes(read_images(mask)) > 0.5).reshape(sigmas[0,...].shape)

    # If there is no mask, make one from the zero set in sigma
    if mask is None:

        mask = np.ones(sigmas[0,...].shape)

    # Loop through m
    for i in np.arange(m):

        # Get sigma for this field
        sigma = sigmas[i,...]

        # Get mask
        mask_current = (sigma != 0) & (~np.isnan(sigma))

        # Take the intersection of the current mask and the previous mask
        mask = np.logical_and(mask, mask_current)
    
    # If tau is none
    if tau is None:

        # If X is none, make it a vector of ones
        if X is None:

            # Make X a list of m (n_sub x 1) matrices with a single entry of 1
            X = [np.ones((n_sub[str(i)],1)) for i in range(m)]

        # If L is none, make it a vector of ones
        if L is None:

            # Make L a list of m (1 x 1) matrices with a single entry of 1
            L = [np.ones((1,1)) for _ in range(m)]

            # Warn the user that L is being set to a vector of ones
            warnings.warn("The contrast vector L was not provided. " \
                          "The analysis will assume we are working with group averages.")
        
        # Check if X and L are lists of numpy arrays
        if not (isinstance(X, list) and isinstance(L, list)):

            # Check if X and L are numpy arrays
            if isinstance(X, np.ndarray) and isinstance(L, np.ndarray):

                # If they are numpy arrays, convert them to lists
                X = [X]
                L = [L]

            # Else raise an error
            else:
                raise TypeError("X and L must be either lists of numpy arrays"\
                                "or numpy arrays themselves.")

        # Make tau a numpy array of zeros with shape (m,1,1)
        tau = np.zeros((m,) + D*(1,))

        # Work out tau by looping through m
        for i in np.arange(m):
            tau_i = np.sqrt(L[i].T @ np.linalg.pinv(X[i].T @ X[i]) @ L[i])[0,0]
            tau[i,...] = tau_i

    # else tau must either be a scalar or numpy array
    else:

        # Check if tau is a scalar or a numpy array
        if isinstance(tau, (int, float, complex)) and not isinstance(tau, bool):

            pass

        elif isinstance(tau, np.ndarray):

            pass
        
        elif isinstance(tau, str):

            # If tau is a filename, read in the file
            tau = cycle_axes(read_images(tau))

        else:
            raise TypeError("tau must be either a scalar or a numpy array")


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

    # If there are no boundary locations, raise an error
    if np.sum(get_bdry_map_combined(cap_muHat, c, mask))==0:

        raise ValueError('No boundary locations found. ' \
                         'Check that the threshold c is ' \
                         'not too high or that the data ' \
                         'is not too noisy.')

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
        for j in np.arange(n_sub[str(i)]):
        
            # Obtain residuals
            resid = cycle_axes(read_images(res_fnames[i][j]))

            # If the mode is sss, we need to standardize the residuals
            if mode == 'sss':

                # Get non-zero sigma values
                sigma = sigmas[i,...].reshape(resid.shape) # MARKER MOVE OUTSIDE LOOP
                sig_mask = (sigma != 0) & (~np.isnan(sigma))

                # Standardize residuals
                resid[sig_mask] = resid[sig_mask]/sigma[sig_mask]

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
    g[mask] = ((muHats[mask]-c)/(tau*sigmas)[mask])

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
        
        # If D is 3, work out image affine
        if D == 3:

            # Read in zeroth resid
            zeroth_resid = nib.load(res_fnames[0][0])

            # Get affine
            affine = zeroth_resid.affine

        # Loop through alpha values saving confidence regions
        for i in np.arange(len(p)):

            # Round p[i] to 2 decimal places
            p[i] = np.round(p[i],2)

            # If D is 3, output nifti image
            if D == 3:

                # Check if previous file exists and if so delete it
                if os.path.exists(os.path.join(out_dir,"Upper_CR_"+str(p[i])+post_fix+".nii")):
                    os.remove(os.path.join(out_dir,"Upper_CR_"+str(p[i])+post_fix+".nii"))
                    warnings.warn(
                        "File Upper_CR_" + str(p[i]) + post_fix + ".nii already exists and will be overwritten."
                    )

                # Save upper confidence region
                addBlockToNifti(os.path.join(out_dir,"Upper_CR_"+str(p[i])+post_fix+".nii"),
                                FcHat_pm_estBdry[i,1,...], np.arange(np.prod(image_dim)),
                                volInd=0,dim=image_dim,aff=affine)

                # Check if previous file exists and if so delete it
                if os.path.exists(os.path.join(out_dir,"Lower_CR_"+str(p[i])+post_fix+".nii")):
                    os.remove(os.path.join(out_dir,"Lower_CR_"+str(p[i])+post_fix+".nii"))
                    warnings.warn(
                        "File Lower_CR_" + str(p[i]) + post_fix + ".nii already exists and will be overwritten."
                    )
                
                # Save filename
                FcHat_plus_files.append(os.path.join(out_dir,"Upper_CR_"+str(p[i])+post_fix+".nii"))
                
                # Save lower confidence region
                addBlockToNifti(os.path.join(out_dir,"Lower_CR_"+str(p[i])+post_fix+".nii"),
                                FcHat_pm_estBdry[i,0,...], np.arange(np.prod(image_dim)),
                                volInd=0,dim=image_dim,aff=affine)
                
                # Save filename
                FcHat_minus_files.append(os.path.join(out_dir,"Lower_CR_"+str(p[i])+post_fix+".nii"))

            # Otherwise output as a numpy array
            else:

                # Check if previous file exists and if so delete it
                if os.path.exists(os.path.join(out_dir,"Upper_CR_"+str(p[i])+post_fix+".npy")):
                    os.remove(os.path.join(out_dir,"Upper_CR_"+str(p[i])+post_fix+".npy"))
                    warnings.warn(
                        "File Upper_CR_" + str(p[i]) + post_fix + ".npy already exists and will be overwritten."
                    )

                # Save upper confidence region
                np.save(os.path.join(out_dir,"Upper_CR_"+str(p[i])+post_fix+".npy"),
                        FcHat_pm_estBdry[i,1,...])
                
                # Save filename
                FcHat_plus_files.append(os.path.join(out_dir,"Upper_CR_"+str(p[i])+post_fix+".npy"))

                # Check if previous file exists and if so delete it
                if os.path.exists(os.path.join(out_dir,"Lower_CR_"+str(p[i])+post_fix+".npy")):
                    os.remove(os.path.join(out_dir,"Lower_CR_"+str(p[i])+post_fix+".npy"))
                    warnings.warn(
                        "File Lower_CR_" + str(p[i]) + post_fix + ".npy already exists and will be overwritten."
                    )

                # Save lower confidence region
                np.save(os.path.join(out_dir,"Lower_CR_"+str(p[i])+post_fix+".npy"),
                        FcHat_pm_estBdry[i,0,...])
                
                # Save filename
                FcHat_minus_files.append(os.path.join(out_dir,"Lower_CR_"+str(p[i])+post_fix+".npy"))
            
        # If D is 3, output nifti image
        if D == 3:
               
            # Check if previous file exists and if so delete it
            if os.path.exists(os.path.join(out_dir,"Estimated_Ac"+post_fix+".nii")):
                os.remove(os.path.join(out_dir,"Estimated_Ac"+post_fix+".nii"))
                warnings.warn(
                    "File Estimated_Ac" + post_fix + ".nii already exists and will be overwritten."
                )
 
            # Save estimated conjunction region
            addBlockToNifti(os.path.join(out_dir,"Estimated_Ac"+post_fix+".nii"),
                            FcHat, np.arange(np.prod(image_dim)),
                            volInd=0,dim=image_dim,aff=affine)
            
            # Save filename
            FcHat_files.append(os.path.join(out_dir,"Estimated_Ac"+post_fix+".nii"))
        
        # Otherwise output as a numpy array
        else:

            # Check if previous file exists and if so delete it
            if os.path.exists(os.path.join(out_dir,"Estimated_Ac"+post_fix+".npy")):
                os.remove(os.path.join(out_dir,"Estimated_Ac"+post_fix+".npy"))
                warnings.warn(
                    "File Estimated_Ac" + post_fix + ".npy already exists and will be overwritten."
                )

            # Save estimated conjunction region
            np.save(os.path.join(out_dir,"Estimated_Ac"+post_fix+".npy"),
                    FcHat)
            
            # Save filename
            FcHat_files.append(os.path.join(out_dir,"Estimated_Ac"+post_fix+".npy"))

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


# Function to generate comparative confidence regions
def generate_conjunction_CRs(mean_fname, sig_fname, res_fnames, out_dir=None, c=None, 
                             p=0.95, mask=None, n_boot=5000, tau=None, X=None, L=None,
                             output=True, mode='sss', method=2, post_fix=None):
    """
    Generate conjunction confidence regions for a set of fields.
    Parameters
    ----------
    mean_fname : str or list of str
        Filename(s) of the mean image(s) for the fields.
    sig_fname : str or list of str
        Filename(s) of the standard deviation image(s) for the fields.
    res_fnames : list of str
        List of filenames for the residuals of the fields.
    out_dir : str, optional
        Output directory for the confidence regions. If None, defaults to the current
        working directory.
    c : float
        Threshold for the conjunction inference. If None, defaults to 0.
    p : float or array_like, optional
        Desired coverage of the confidence regions. If a single float, it is used for
        all fields. If an array, it should have the same length as the number of fields.
        Defaults to 0.95.
    mask : array_like, optional 
        Mask for the data. If None, defaults to a mask of ones with the same shape as
        the mean image.
    n_boot : int, optional
        Number of bootstrap samples to use for estimating the confidence regions.
        Defaults to 5000.
    tau : float or array_like, optional
        Image or scalar of tau_n values. If None, defaults to sqrt(n).
    X : array_like, optional
        Design matrix in a regression context. If None, a vector of ones is used,
        and a simple mean model is assumed.
    L : array_like, optional
        Contrast matrix in a regression model. If None, a (1,1) matrix with a single
        entry of 1 is used, and a simple mean model is assumed.
    output : bool, optional
        Whether to output the confidence regions as nifti files. If False, the
        confidence regions are returned as numpy arrays. Defaults to True.
    mode : str, optional
        Type of confidence regions to generate. Can be 'sss' for the method of Sommerfeld,
        Sain and Schwartzman (2018), or 'cohens' for the cohens d method of Bowring et al.
        (2019). Defaults to 'sss'.
    method : int, optional
        Method to use for the cohens d method. Can be 1, 2 or 3. The three algorithms
        are those listed in the following manuscript:
            https://doi.org/10.1016/j.neuroimage.2020.117477
        (c.f. Section 2.6). Method 2 is used by default following the recommendations of
        the above reference. This argument is not needed for 'sss' mode.
    post_fix : str, optional
        Post-fix to add to the output files. If None, no post-fix is added.
        Defaults to None.
    Returns
    -------
    FcHat_plus_files : list of str
        List of filenames for the upper confidence regions.
    FcHat_minus_files : list of str
        List of filenames for the lower confidence regions.
    FcHat_files : list of str
        List of filenames for the estimated conjunction region.
    a_estBdry : array_like
        Estimated boundary values for the conjunction region.
    """
    
    # Call to generate_CRs to generate conjunction confidence regions
    FcHat_plus_files, FcHat_minus_files, FcHat_files, a_estBdry = generate_CRs(
        mean_fname, sig_fname, res_fnames, out_dir=out_dir, c=c, p=p, mask=mask,
        n_boot=n_boot, tau=tau, X=X, L=L, output=output, mode=mode, method=method,
        post_fix=post_fix
    )

    # Return results
    return FcHat_plus_files, FcHat_minus_files, FcHat_files, a_estBdry


# Function to generate comparative confidence regions
def generate_set_diff_CRs(mean_fname, sig_fname, res_fnames, out_dir=None, c=None, 
                          p=0.95, mask=None, n_boot=5000, tau=None, X=None, L=None,
                          output=True, mode='sss', method=2, post_fix=None):
    
    """
    Generate comparative confidence regions for a set of fields.
    Parameters
    ----------
    mean_fname : list of str
        List of filenames for the mean images of the two fields.
    sig_fname : str or list of str
        Filename(s) of the standard deviation image(s) for the fields.
    res_fnames : list of str
        List of filenames for the residuals of the fields.
    out_dir : str, optional
        Output directory for the confidence regions. If None, defaults to the current
        working directory.
    c : float
        Threshold for the comparative inference. If None, defaults to 0.
    p : float or array_like, optional
        Desired coverage of the confidence regions. If a single float, it is used for
        both fields. Defaults to 0.95.
    mask : array_like, optional
        Mask for the data. If None, defaults to a mask of ones with the same shape as
        the mean image.
    n_boot : int, optional
        Number of bootstrap samples to use for estimating the confidence regions.
        Defaults to 5000.
    tau : float or array_like, optional
        Image or scalar of tau_n values. If None, defaults to sqrt(n).
    X : array_like, optional
        Design matrix in a regression context. If None, a vector of ones is used,
        and a simple mean model is assumed.
    L : array_like, optional
        Contrast matrix in a regression model. If None, a matrix with a first entry
        entry of 1 is used, and a simple mean model is assumed.
    output : bool, optional
        Whether to output the confidence regions as nifti files. If False, the
        confidence regions are returned as numpy arrays. Defaults to True.
    mode : str, optional
        Type of confidence regions to generate. Can be 'sss' for the method of Sommerfeld,
        Sain and Schwartzman (2018), or 'cohens' for the cohens d method of Bowring et al.
        (2019). Defaults to 'sss'.
    method : int, optional
        Method to use for the cohens d method. Can be 1, 2 or 3. The three algorithms
        are those listed in the following manuscript:
            https://doi.org/10.1016/j.neuroimage.2020.117477
        (c.f. Section 2.6). Method 2 is used by default following the recommendations of
        the above reference. This argument is not needed for 'sss' mode.
    post_fix : str, optional
        Post-fix to add to the output files. If None, no post-fix is added.
        Defaults to None.
    Returns
    -------
    FcHat_plus_files : list of str
        List of filenames for the upper confidence regions.
    FcHat_minus_files : list of str
        List of filenames for the lower confidence regions.
    FcHat_files : list of str   
        List of filenames for the estimated comparative region.
    a_estBdry : array_like  
        Estimated boundary values for the comparative region.
    """

    
    # Check if mean_fname is a list
    if isinstance(mean_fname, list):

        # Set m equal to the length of mean_fname
        m = len(mean_fname)

        # If m is not equal to 2, raise an error
        if m != 2:
            raise ValueError('Comparative confidence regions can only be generated for two fields.')

    # If not error
    else:
        raise TypeError('mean_fname must be a list of filenames.')
    
    # Read in second mean volume
    muHat2 = read_image(mean_fname[1])

    # Get D
    D = len(muHat2.shape)

    # Compute c minus muHat2
    muHat2 = 2*c - muHat2

    # Save muHat2
    if D == 3:

        # Create modified filename
        new_mean_fname = os.path.join(mean_fname[1].replace('.nii', '_set_diff.nii'))
        
        # Must save as nifti file
        addBlockToNifti(new_mean_fname, muHat2, 
                        np.arange(np.prod(muHat2.shape)), volInd=0, dim=muHat2.shape)
        
        # Set mean_fname to new filename
        mean_fname[1] = new_mean_fname
        
    # Otherwise save as numpy array
    else:

        # Create modified filename
        new_mean_fname = os.path.join(mean_fname[1].replace('.npy', '_set_diff.npy'))

        # Save as numpy array
        np.save(new_mean_fname, muHat2)

        # Set mean_fname to new filename
        mean_fname[1] = new_mean_fname

    # Call generate_CRs with modified mean_fname
    FcHat_plus_files, FcHat_minus_files, FcHat_files, a_estBdry = generate_CRs(
        mean_fname, sig_fname, res_fnames, out_dir=out_dir, c=c, p=p, mask=mask,
        n_boot=n_boot, tau=tau, X=X, L=L, output=output, mode=mode, method=method,
        post_fix=post_fix
    )

    # Remove new mean file
    if os.path.exists(new_mean_fname):
        os.remove(new_mean_fname)

    # Return results
    return FcHat_plus_files, FcHat_minus_files, FcHat_files, a_estBdry


    