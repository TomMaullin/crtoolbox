import os
import numpy as np
import nibabel as nib
from crtoolbox.lib.fileio import read_image, remove_files
from crtoolbox.tests.generate_ni_data import addBlockToNifti
from crtoolbox.lib.regression import *


"""
Transform a list of data images into cohen's d residual images.

Parameters:
-----------
    data_fnames : list
        list of filenames for data
    X : array
        design matrix
    out_dir : str
        output directory
    method: str
        Number between 1 and 3 indicating which method to use for calculating
        cohen's d residuals. The three algoirthms are those listed in the 
        following manuscript:
            https://doi.org/10.1016/j.neuroimage.2020.117477
        (c.f. Section 2.6)
        Method 2 is used by default following the recommendations of the above
        reference.

Returns:
--------
    d_fname : str
        filename for cohen's d image
    cohen_fnames : list
        list of filenames for cohen's d residual images
    cohen_sigma_fname : str
        filename for cohen's d residual standard deviation image (only returned
        if method=2)
"""
def cohens(data_fnames, X, out_dir, method=2):

    # Get mean and standard deviation images
    _, mean_fname, std_fname, resid_files, _ = regression(data_fnames, X, out_dir)  # MARKER - Currently not accounting for choice of contrast

    # Remove the residual files
    remove_files(resid_files)  

    # Check if output directory exists
    if not os.path.exists(out_dir):
            
        # Create output directory
        os.mkdir(out_dir)

    # Work out number of data images
    n = len(data_fnames)

    # Read in mean image
    mean_img = read_image(mean_fname)

    # Read in standard deviation image
    std_img = read_image(std_fname)

    # Read in a single image to get image size
    img = read_image(data_fnames[0])

    # Get image size
    img_size = img.shape

    # Make sure the mean and standard deviation images are the same size
    mean_img = mean_img.reshape(img_size)
    std_img = std_img.reshape(img_size)

    # Get the dimension of the image
    D = len(img_size)

    # Initialize list of cohen's d residual filenames
    cohen_fnames = []

    # Get image of where sigma is non-zero
    std_nonzero = std_img > 0

    # Get non-zero values of mean image and standard deviation image
    mean_img_nonzero = mean_img[std_nonzero]
    std_img_nonzero = std_img[std_nonzero]

    # Compute cohen's d
    cohens_d = np.zeros(img_size)
    cohens_d[std_nonzero] = mean_img_nonzero / std_img_nonzero

    # Save cohen's d to file
    if D == 3:

        # Write cohen's d to file
        addBlockToNifti(os.path.join(out_dir,"cohen_d.nii"), cohens_d,
                        np.arange(np.prod(img_size)), volInd=0,dim=img_size)
        
        # Add filename to list
        d_fname = os.path.join(out_dir,"cohen_d.nii")

    else:

        # Write cohen's d to file
        np.save(os.path.join(out_dir,"cohen_d.npy"), cohens_d)

        # Add filename to list
        d_fname = os.path.join(out_dir,"cohen_d.npy")

    # Check if we are using method 1
    if method == 1:

        # Loop through data images
        for i in range(0, n):

            # Read in data image
            data_img = read_image(data_fnames[i])

            # Get epsilon
            epsilon = data_img[std_nonzero] - mean_img_nonzero

            # Get cohen's d residual
            resid = epsilon / std_img_nonzero - (mean_img_nonzero / (2*std_img_nonzero)*((epsilon / std_img_nonzero)**2 - 1))

            # Normalise the residual
            resid_tilde = np.zeros(img_size)
            resid_tilde[std_nonzero] = resid / np.sqrt(1 + cohens_d[std_nonzero]**2/2)

            # Output residual to file
            if D == 3:

                # Write residual to file
                addBlockToNifti(os.path.join(out_dir,"cohen_" + str(i) + ".nii"), resid_tilde,
                                np.arange(np.prod(img_size)), volInd=0,dim=img_size)
                
                # Add filename to list
                cohen_fnames.append(os.path.join(out_dir,"cohen_" + str(i) + ".nii"))

            else:

                # Write residual to file
                np.save(os.path.join(out_dir,"cohen_" + str(i) + ".npy"), resid_tilde)

                # Add filename to list
                cohen_fnames.append(os.path.join(out_dir,"cohen_" + str(i) + ".npy"))

        # Set cohen sigma filename to none
        cohen_sigma_fname = None

    # Check if we are using method 2
    elif method == 2:

        # Loop through data images to compute mean of residuals
        for i in range(0, n):

            # Read in data image
            data_img = read_image(data_fnames[i])

            # Get epsilon
            epsilon = data_img[std_nonzero] - mean_img_nonzero

            # Get cohen's d residual
            resid = epsilon / std_img_nonzero - (mean_img_nonzero / (2*std_img_nonzero)*((epsilon / std_img_nonzero)**2 - 1))

            # Add to running sum
            if i == 0:

                # Initialize running sum
                resid_sum = resid

                # Initialize running sum of squares
                resid_sum_sq = resid**2

            else:

                # Add to running sum
                resid_sum += resid

                # Add to running sum of squares
                resid_sum_sq += resid**2

        # Get mean of residuals
        resid_mean = resid_sum / n

        # Get standard deviation of residuals
        resid_sigma = np.sqrt((resid_sum_sq - n * resid_mean**2) / (n - 1))

        # Output resid sigma to file
        resid_sigma_vol = np.zeros(img_size)
        resid_sigma_vol[std_nonzero] = resid_sigma

        # Output resid sigma to file
        if D == 3:

            # Write residual to file
            addBlockToNifti(os.path.join(out_dir,"cohen_sigma.nii"), resid_sigma_vol,
                            np.arange(np.prod(img_size)), volInd=0,dim=img_size)
            
            # Save filename
            cohen_sigma_fname = os.path.join(out_dir,"cohen_sigma.nii")

        else:

            # Write residual to file
            np.save(os.path.join(out_dir,"cohen_sigma.npy"), resid_sigma_vol)
            
            # Save filename
            cohen_sigma_fname = os.path.join(out_dir,"cohen_sigma.npy")

        # Loop through data images to compute normalised residuals
        for i in range(0, n):

            # Read in data image
            data_img = read_image(data_fnames[i])

            # Get epsilon
            epsilon = data_img[std_nonzero] - mean_img_nonzero

            # Get cohen's d residual
            resid = epsilon / std_img_nonzero - (mean_img_nonzero / (2*std_img_nonzero)*((epsilon / std_img_nonzero)**2 - 1))

            # Get the normalized residual
            resid_tilde = np.zeros(img_size)
            resid_tilde[std_nonzero] = resid / resid_sigma

            # Output residual to file
            if D == 3:

                # Write residual to file
                addBlockToNifti(os.path.join(out_dir,"cohen_" + str(i) + ".nii"), resid_tilde,
                                np.arange(np.prod(img_size)), volInd=0,dim=img_size)
                
                # Add filename to list
                cohen_fnames.append(os.path.join(out_dir,"cohen_" + str(i) + ".nii"))

            else:

                # Write residual to file
                np.save(os.path.join(out_dir,"cohen_" + str(i) + ".npy"), resid_tilde)

                # Add filename to list
                cohen_fnames.append(os.path.join(out_dir,"cohen_" + str(i) + ".npy"))

    # Check if we are using method 3
    elif method == 3:

        # Compute constants
        alpha_star = 1/np.sqrt((n * (8 * n**2 - 17 * n + 11)) / ((n - 3) * (4 * n - 5)**2))
        beta_star = np.sqrt((n * (n - 1)) / (n - 3))/np.sqrt((8 * n**2 - 17 * n + 11) / ((n - 3) * (4 * n - 5)**2))
        # b_star = np.sqrt((n * (8 * n**2 - 17 * n + 11)) / ((n - 3) * (4 * n - 5)**2))
    

        # Loop through data images
        for i in range(0, n):

            # Read in data image
            data_img = read_image(data_fnames[i])

            # Get epsilon
            epsilon = data_img[std_nonzero] - mean_img_nonzero

            # Get variance stabilized cohen's d residual
            resid_tilde = np.zeros(img_size)
            resid_tilde[std_nonzero] = (epsilon / std_img_nonzero - (mean_img_nonzero / (2*std_img_nonzero)*((epsilon / std_img_nonzero)**2 - 1)))*alpha_star*beta_star/np.sqrt(1 + (beta_star*epsilon)**2/2)

            # Output residual to file
            if D == 3:

                # Write residual to file
                addBlockToNifti(os.path.join(out_dir,"cohen_" + str(i) + ".nii"), resid_tilde,
                                np.arange(np.prod(img_size)), volInd=0,dim=img_size)
                
                # Add filename to list
                cohen_fnames.append(os.path.join(out_dir,"cohen_" + str(i) + ".nii"))

            else:

                # Write residual to file
                np.save(os.path.join(out_dir,"cohen_" + str(i) + ".npy"), resid_tilde)

                # Add filename to list
                cohen_fnames.append(os.path.join(out_dir,"cohen_" + str(i) + ".npy"))

        # Set cohen sigma filename to None
        cohen_sigma_fname = None

    # Otherwise throw error
    else:
        
        # Throw error
        raise ValueError("Invalid method specified. Must be 1, 2 or 3.")   

    # We no longer need the original mean and std files
    remove_files(mean_fname)
    remove_files(std_fname) 

    # Return list of cohen's d residual filenames
    return d_fname, cohen_fnames, cohen_sigma_fname






