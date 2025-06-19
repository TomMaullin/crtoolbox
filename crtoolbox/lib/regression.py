import os
import numpy as np
import nibabel as nib
from crtoolbox.lib.fileio import read_image
from crtoolbox.tests.generate_ni_data import addBlockToNifti

""" 
Read images from a list of filenames one by one and compute their mean and standard
deviation. If the parameter mean_zero is set to True then we assume that the images
are already mean zero.

Parameters:
-----------
fnames : list
    list of filenames
mean_zero: bool
    If true then we assume the images are already mean zero.

Returns:
--------
mean_img : array
    mean image
std_img : array
    standard deviation image
"""
def read_images_mean_std(fnames, mean_zero=False):

    # If we are not assuming the images are already mean zero
    if not mean_zero:

        # Loop through files
        for i in range(0, len(fnames)):

            # Read in image
            img = read_image(fnames[i])

            # Check if we are on the first image
            if i == 0:
                    
                # Initialize mean array
                sum_img = np.array(img, dtype=np.float64)

                # Record image size
                img_size = img.shape

            else:

                # Take sum
                sum_img = sum_img + img

                # Check if image is the same size as the first image
                if img.shape != img_size:

                    # Raise error
                    raise ValueError("Error: Images are not all the same size")

                    # Exit
                    exit()


        # Compute mean
        mean_img = sum_img / len(fnames)

        # Loop through files
        for i in range(0, len(fnames)):

            # Read in image
            img = read_image(fnames[i])

            # Check if we are on the first image
            if i == 0:

                # Initialize array
                sum_sq_img = np.array((img - mean_img)**2, dtype=np.float64)

            else:

                # Take sum
                sum_sq_img = sum_sq_img + (img - mean_img)**2

        # Compute standard deviation
        std_img = np.sqrt(sum_sq_img / len(fnames))

    else:

        # Loop through files
        for i in range(0, len(fnames)):

            # Read in image
            img = read_image(fnames[i])

            # Check if we are on the first image
            if i == 0:

                # Initialize array
                sum_sq_img = np.array(img**2, dtype=np.float64)

            else:

                # Take sum
                sum_sq_img = sum_sq_img + img**2

        # Compute standard deviation
        std_img = np.sqrt(sum_sq_img / len(fnames))

        # Compute mean
        mean_img = np.zeros(std_img.shape)

    return mean_img, std_img


"""
Read images from a list of filenames in chunks and perform
regression on each chunk. The regression is performed using
the following equation:

    y = X * beta + epsilon

where y is the image, X is the design matrix, beta is the 
regression coefficients, and epsilon is the error. The
regression coefficients are computed using the following
equation:

    beta = (X'X)^-1 X'y

where X' is the transpose of X and ^-1 is the inverse.

Parameters:
-----------
yfiles : list
    list of filenames for the images to be regressed
X : array
    design matrix
out_dir : str
    directory to save the regression results
L: array
    contrast vector for the regression (default has 1 
    for the first row and 0 for the rest). Input is
    assumed to be a row vector of shape (1, p).
chunk_size : int
    number of images to read in at a time
post_fix : str
    post fix to add to the output filenames (default is 
    empty string)

Returns:
--------
files : list
    list of filenames for the regression results. The filenames
    are in the order of: beta, sigma, and epsilon.
"""
def regression(yfiles, X, out_dir, L=None, chunk_size=20, post_fix=""):

    # Check if output directory exists
    if not os.path.exists(out_dir):
        # Create output directory
        os.makedirs(out_dir)

    # Check if yfiles is a list
    if not isinstance(yfiles, list):
        # Raise error
        raise ValueError("yfiles must be a list of filenames")
    
    # Check if X is a numpy array
    if not isinstance(X, np.ndarray):
        # Raise error
        raise ValueError("X must be a numpy array")
    
    # Check if X is 2D
    if len(X.shape) != 2:
        # Raise error
        raise ValueError("X must be a 2D numpy array")
    
    # If the post fix is not empty
    if post_fix != "":
        # Add an underscore to the post fix
        post_fix = "_" + post_fix
    
    # Check if L is None
    if L is None:
        # Create L as a vector with 1 for the first column 
        # and 0 for the rest
        L = np.array([[1] + [0]*(X.shape[-1]-1)])
    # If L is 1D make it a row vector
    elif L.ndim == 1:
        # Convert L to a row vector
        L = np.array(L).reshape(1, -1)
    # If L is 2D but a column vector
    elif L.shape[1] == 1:
        # Transpose L to ensure it is a row vector
        L = np.array(L).T

    # Get number of images
    n_imgs = len(yfiles)

    # Get number of rows in design matrix
    n = X.shape[0]

    # Check if n_imgs matches number of rows in design matrix
    if n_imgs != n:
            
        # Raise error
        raise ValueError("Number of images does not match number of rows in design matrix")

    # Get number of columns in design matrix
    p = X.shape[1]

    # Read in a single image to get image size
    img = read_image(yfiles[0])

    # Get image size
    img_size = img.shape

    # Get the dimension of the image
    D = len(img_size)

    # Create an array to store na sum
    na_sum = np.zeros(img_size)

    # Read in images in chunks
    for i in range(0, n_imgs, chunk_size):

        # Get end index
        end = min(i + chunk_size, n_imgs)

        # Get number of images in chunk
        n_chunk = end - i

        # Initialize arrays
        y_chunk = np.zeros(img_size + (n_chunk, 1))
        X_chunk = np.zeros(tuple((1 for i in range(D))) + (n_chunk, p))

        # Loop through images in chunk        
        for j in range(0, n_chunk):

            # Read in image
            img = read_image(yfiles[i + j])

            # Add image to array
            y_chunk[..., j, 0] = img

            # Add design matrix to array
            X_chunk[..., j, :] = X[i + j, :]

        # Compute X'X for chunk
        XtX_chunk = X_chunk.swapaxes(-1,-2) @ X_chunk

        # Compute X'y for chunk
        Xty_chunk = X_chunk.swapaxes(-1,-2) @ y_chunk

        # Handle NaNs in Xty_chunk
        Xty_chunk[np.isnan(Xty_chunk)] = 0  # Handle NaNs

        # Compute y'y for chunk
        yty_chunk = y_chunk.swapaxes(-1,-2) @ y_chunk

        # Handle NaNs in yty_chunk
        yty_chunk[np.isnan(yty_chunk)] = 0  # Handle NaNs

        # Compute na_sum for chunk
        na_sum_chunk = np.squeeze(np.sum(np.isnan(y_chunk), axis=-2, keepdims=True))

        # Add X'X, X'y and y'y for chunk to total
        if i == 0:

            # Initialize arrays
            XtX = np.zeros(XtX_chunk.shape)
            Xty = np.zeros(Xty_chunk.shape)
            yty = np.zeros(yty_chunk.shape)

        # Add to total
        XtX = XtX + XtX_chunk
        Xty = Xty + Xty_chunk
        yty = yty + yty_chunk

        # Add to na_sum
        na_sum = na_sum + na_sum_chunk


    # Missing data threshold
    mt = 0.05

    # Create mask using na_sum
    mask = (na_sum < n_imgs*mt)

    # Compute beta
    beta = np.linalg.pinv(XtX) @ Xty

    # Apply mask to beta
    beta[~mask] = np.nan

    # Compute sum of squared errors
    ete = yty - beta.swapaxes(-1,-2) @ Xty

    # Compute sigma
    sigma = np.sqrt(ete / n)

    # Apply mask to sigma
    sigma[~mask] = np.nan

    # Empty stores for filenames
    beta_files = []
    resid_files = []
    sigma_file = []

    # Pad contrast vector dimensions
    if len(L.shape) < len(beta.shape):

        # Make sure L has extra dimensions for broadcasting
        L = L.reshape((1,) * (len(beta.shape) - len(L.shape)) + L.shape)

    # Compute beta contrast
    contrast = L @ beta

    # If the image is 3D then output nifti image
    if D == 3:

        # Write sigma to file
        addBlockToNifti(os.path.join(out_dir,"sigma" + post_fix + ".nii"), sigma,
                        np.arange(np.prod(img_size)), volInd=0,dim=img_size)
        
        # Add filename to list
        sigma_file = os.path.join(out_dir,"sigma" + post_fix + ".nii")

        # Write mask to file
        addBlockToNifti(os.path.join(out_dir,"mask" + post_fix + ".nii"), np.abs(sigma) > 1e-8,
                        np.arange(np.prod(img_size)), volInd=0,dim=img_size)
        
        # Save mask filename
        mask_file = os.path.join(out_dir,"mask" + post_fix + ".nii")
        
        # Write na_sum to file
        # addBlockToNifti(os.path.join(out_dir,"na_sum" + post_fix + ".nii"), na_sum,
        #                 np.arange(np.prod(img_size)), volInd=0,dim=img_size)
        
        # Write contrast to file
        addBlockToNifti(os.path.join(out_dir,"contrast" + post_fix + ".nii"), contrast[..., 0],
                        np.arange(np.prod(img_size)), volInd=0,dim=img_size)
        
        # Save contrast filename
        contrast_file = os.path.join(out_dir,"contrast" + post_fix + ".nii")

    # Otherwise output as a numpy array
    else:

        # Write sigma to file
        np.save(os.path.join(out_dir,"sigma" + post_fix + ".npy"), sigma)

        # Add filename to list
        sigma_file = os.path.join(out_dir,"sigma" + post_fix + ".npy")

        # Write sigma to file
        np.save(os.path.join(out_dir,"mask" + post_fix + ".npy"), sigma != 0)

        # Save mask filename
        mask_file = os.path.join(out_dir,"mask" + post_fix + ".npy")

        # Write na_sum to file
        # np.save(os.path.join(out_dir,"na_sum" + post_fix + ".npy"), na_sum)

        # Write contrast to file
        np.save(os.path.join(out_dir,"contrast" + post_fix + ".npy"), contrast[..., 0])

        # Save contrast filename
        contrast_file = os.path.join(out_dir,"contrast" + post_fix + ".npy")

    # Loop through beta coefficients
    for i in range(0, p):

        # If the image is 3D then output nifti image
        if D == 3:

            # Write beta coefficient to file
            addBlockToNifti(os.path.join(out_dir,"betahat"+str(i)+ post_fix + ".nii"), 
                            beta[..., i, 0], np.arange(np.prod(img_size)), 
                            volInd=0,dim=img_size)
            
            # Add filenames to list
            beta_files.append(os.path.join(out_dir,"betahat"+str(i) + post_fix + ".nii"))
            
        # Otherwise output as a numpy array
        else:

            # Write beta coefficient to file
            np.save(os.path.join(out_dir,"betahat"+str(i) + post_fix + ".npy"), beta[..., i, 0])

            # Add filenames to list
            beta_files.append(os.path.join(out_dir,"betahat"+str(i) + post_fix + ".npy"))

    # Loop through images creating residuals
    for i in range(0, n_imgs):

        # Read in image
        img = read_image(yfiles[i])

        # Compute residuals
        res = (img - (X[..., i:(i+1), :] @ beta)[..., 0, 0])

        # Apply mask to residuals
        res[~mask] = np.nan

        # If the image is 3D then output nifti image
        if D == 3:

            # Write residuals to file
            addBlockToNifti(os.path.join(out_dir,"res"+str(i) + post_fix + ".nii"), 
                            res, np.arange(np.prod(img_size)), 
                            volInd=0,dim=img_size)
            
            # Add filenames to list
            resid_files.append(os.path.join(out_dir,"res"+str(i) + post_fix + ".nii"))

        # Otherwise output as a numpy array
        else:

            # Write residuals to file
            np.save(os.path.join(out_dir,"res"+str(i) + post_fix + ".npy"), res)

            # Add filenames to list
            resid_files.append(os.path.join(out_dir,"res"+str(i) + post_fix + ".npy"))

    # Check if there is only one beta file 
    if len(beta_files) == 1:

        # Change beta_files to a string
        beta_files = beta_files[0]
        
    # Return filenames
    return contrast_file, beta_files, sigma_file, resid_files, mask_file

        