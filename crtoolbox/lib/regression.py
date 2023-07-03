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
chunk_size : int
    number of images to read in at a time

Returns:
--------
files : list
    list of filenames for the regression results. The filenames
    are in the order of: beta, sigma, and epsilon.
"""
def regression(yfiles, X, out_dir, chunk_size=20):

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

        # Compute y'y for chunk
        yty_chunk = y_chunk.swapaxes(-1,-2) @ y_chunk

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

    # Compute beta
    beta = np.linalg.pinv(XtX) @ Xty

    # Compute sum of squared errors
    ete = yty - beta.swapaxes(-1,-2) @ Xty

    # Compute sigma
    sigma = np.sqrt(ete / n)

    # Empty stores for filenames
    beta_files = []
    resid_files = []
    sigma_file = []

    # If the image is 3D then output nifti image
    if D == 3:

        # Write sigma to file
        addBlockToNifti(os.path.join(out_dir,"sigma.nii"), sigma,
                        np.arange(np.prod(img_size)), volInd=0,dim=img_size)
        
        # Add filename to list
        sigma_file = os.path.join(out_dir,"sigma.nii")

    # Otherwise output as a numpy array
    else:

        # Write sigma to file
        np.save(os.path.join(out_dir,"sigma.npy"), sigma)

        # Add filename to list
        sigma_file = os.path.join(out_dir,"sigma.npy")

    # Loop through beta coefficients
    for i in range(0, p):

        # If the image is 3D then output nifti image
        if D == 3:

            # Write beta coefficient to file
            addBlockToNifti(os.path.join(out_dir,"betahat"+str(i)+".nii"), 
                            beta[..., i, 0], np.arange(np.prod(img_size)), 
                            volInd=0,dim=img_size)
            
            # Add filenames to list
            beta_files.append(os.path.join(out_dir,"betahat"+str(i)+".nii"))
            
        # Otherwise output as a numpy array
        else:

            # Write beta coefficient to file
            np.save(os.path.join(out_dir,"betahat"+str(i)+".npy"), beta[..., i, 0])

            # Add filenames to list
            beta_files.append(os.path.join(out_dir,"betahat"+str(i)+".npy"))

    # Loop through images creating residuals
    for i in range(0, n_imgs):

        # Read in image
        img = read_image(yfiles[i])

        # Compute residuals
        res = (img - (X[..., i:(i+1), :] @ beta)[..., 0, 0])

        # If the image is 3D then output nifti image
        if D == 3:

            # Write residuals to file
            addBlockToNifti(os.path.join(out_dir,"res"+str(i)+".nii"), 
                            res, np.arange(np.prod(img_size)), 
                            volInd=0,dim=img_size)
            
            # Add filenames to list
            resid_files.append(os.path.join(out_dir,"res"+str(i)+".nii"))

        # Otherwise output as a numpy array
        else:

            # Write residuals to file
            np.save(os.path.join(out_dir,"res"+str(i)+".npy"), res)

            # Add filenames to list
            resid_files.append(os.path.join(out_dir,"res"+str(i)+".npy"))

    # Check if there is only one beta file 
    if len(beta_files) == 1:

        # Change beta_files to a string
        beta_files = beta_files[0]
        
    # Return filenames
    return beta_files, sigma_file, resid_files  
        