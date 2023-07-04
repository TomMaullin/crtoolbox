import warnings as w
# This warning is caused by numpy updates and should
# be ignored for now.
w.simplefilter(action = 'ignore', category = FutureWarning)
import numpy as np
import nibabel as nib
import sys
import os
import glob
import shutil
import yaml
from crtoolbox.lib.fileio import *
import time
import pandas as pd
from scipy import ndimage

def generate_data(n, p, OutDir, dim=np.array([100,100,100]), mask_type='fixed', sigma=3):
    """
    Generates simulated data with the specified dimensions and other parameters.
    
    Parameters
    ----------
    n : int
        Number of observations.
    p : int
        Number of parameters.
    dim : numpy array
        Dimensions of data to be generated. Must be given as a numpy array.
    OutDir : str
        Output directory path where the generated data will be saved.
    mask_type : str
        Type of mask to be used. Must be either "random" or "fixed".
    sigma : float
        Standard deviation of the noise to be added to the data.

    Returns
    -------
    None
        This function does not return any value. It saves the generated data in the specified output directory.
    """

    # Check if data directory exists
    if not os.path.exists(OutDir):
        os.mkdir(OutDir)

    # Make sure in numpy format (added 20 for smoothing)
    origdim = np.array(dim)
    dim = origdim + 20

    # -------------------------------------------------
    # Design parameters
    # -------------------------------------------------

    # fwhm for smoothing
    fwhm = 5

    # Number of voxels
    v = np.prod(dim)

    # -------------------------------------------------
    # Obtain design matrix
    # -------------------------------------------------

    # Fixed effects design matrix
    X = get_X(n, p)
    
    # -------------------------------------------------
    # Obtain beta parameter vector
    # -------------------------------------------------

    # Get beta
    beta = get_beta(p, dim)

    # -----------------------------------------------------
    # Obtain Y
    # -----------------------------------------------------

    # Loop through subjects generating nifti images
    for i in np.arange(n):    

        # Initialize Yi to Xi times beta
        Yi = (X[...,i:(i+1),:] @ beta)[...,0,0]

        # Get epsiloni
        epsiloni = get_epsilon(v, 1, sigma).reshape(dim)

        # Smooth epsiloni
        epsiloni = smooth_data(epsiloni, 3, [fwhm]*3, trunc=6, scaling='kernel').reshape(dim)

        # Add epsilon to Yi
        Yi = Yi + epsiloni

        # Obtain mask
        if mask_type == 'random':
            # Get random mask
            mask = get_random_mask(dim).reshape(Yi.shape)
        elif mask_type == 'fixed':
            # Get mask from file
            mask = nib.load(os.path.join(os.path.dirname(__file__),'mask.nii')).get_data()
        else:
            raise ValueError('mask_type must be either "random" or "fixed".')

        # Mask Yi
        Yi = Yi*mask

        # Truncate off (handles smoothing edge effects)
        Yi = Yi[10:(dim[0]-10),10:(dim[1]-10),10:(dim[2]-10)]

        # Output Yi
        addBlockToNifti(os.path.join(OutDir,"Y"+str(i)+".nii"), Yi, np.arange(np.prod(origdim)), volInd=0,dim=origdim)

    # Reshape X
    X = X.reshape(n,p)

    # -----------------------------------------------------
    # Save beta
    # -----------------------------------------------------

    # Mask beta
    beta = beta*mask.reshape(mask.shape + (1,1))

    # Truncate beta
    beta = beta[10:(dim[0]-10),10:(dim[1]-10),10:(dim[2]-10),:,:]

    # Loop through beta values
    for i in np.arange(p):

        # Save beta to nifti
        addBlockToNifti(os.path.join(OutDir,"beta"+str(i)+".nii"), beta[...,i,0], np.arange(np.prod(origdim)), volInd=0,dim=origdim)

    # -----------------------------------------------------
    # Yfiles storage
    # -----------------------------------------------------

    # Create empty list to store Y filenames
    Yfiles = []

    # Loop through listing y files in text file
    for i in np.arange(n):

        # Add current Y filename to list
        Yfiles.append(os.path.join(OutDir,"Y"+str(i)+".nii"))

    # -----------------------------------------------------
    # betafiles storage
    # -----------------------------------------------------

    # Create empty list to store beta filenames
    betafiles = []

    # Loop through listing beta files in text file
    for i in np.arange(p):

        # Add current beta filename to list
        betafiles.append(os.path.join(OutDir,"beta"+str(i)+".nii"))


    return Yfiles, betafiles, X


def get_random_mask(dim):
    """
    Generates a random mask with the specified dimensions.
    
    Parameters
    ----------
    dim: numpy array
        Dimensions of data to be generated. Must be given as a numpy array.

    Returns
    -------
    mask: numpy array
        Random mask with the specified dimensions.
    """

    # FWHM
    fwhm = 10

    # Load analysis mask
    mask = nib.load(os.path.join(os.path.dirname(__file__),'mask.nii')).get_data()

    # Add some noise and smooth
    mask = smooth_data(mask + 8*np.random.randn(*(mask.shape)), 3, [fwhm]*3)

    # Re-threshold (this has induced a bit of randomness in the mask shape)
    mask = 1*(mask > 0.6)

    return(mask)


def get_X(n,p):
    """
    Generates a random design matrix with the specified dimensions.

    Parameters
    ----------
    n : int
        Number of observations.
    p : int
        Number of parameters.

    Returns
    -------
    X : numpy array for design matrix
    """

    # Generate random X.
    X = np.random.uniform(low=-0.5,high=0.5,size=(n,p))
    
    # Make the first column an intercept
    X[:,0]=1

    # Reshape to dimensions for broadcasting
    X = X.reshape(1, n, p)

    # Return X
    return(X)


def get_beta(p, dim):
    """
    Generates a random beta parameter vector with the specified dimensions.

    Parameters
    ----------
    p : int
        Number of parameters.
    dim : numpy array
        Dimensions of data to be generated. Must be given as a numpy array.

    Returns
    -------
    beta : numpy array for beta parameter vector
    """

    # fwhm
    fwhm = 10

    # Loop through beta values
    for i in np.arange(p):

        # Without assuming the number of dimension, for each dimension,
        # generate a random coordinate near the center
        # if i == 0:
        #     center = np.array([50,50,70])
        # else:
        center = np.array([np.random.uniform(low=0.4*dim[i],high=0.6*dim[i]) for i in np.arange(len(dim))])

        # Generate a random radius between 5 and 12
        # if i == 0:
        #     r = 10
        # else:
        r = np.random.uniform(low=5,high=10)

        # Create an ogrid
        ogrid = np.meshgrid(*[np.arange(d) for d in dim], indexing='ij')

        # Calculate distance from the center
        distance = sum((g - c)**2 for g, c in zip(ogrid, center))

        # Create the sphere and rescale
        beta_p = np.array(np.sqrt(distance) < r, dtype='float')

        # Smooth the sphere
        beta_p = smooth_data(beta_p, 3, [fwhm]*3, trunc=6, scaling='max')*3
    
        # Concatenate betas along additional first axis
        if i == 0:
            beta = beta_p.reshape(*beta_p.shape, 1, 1)
        else:
            beta = np.concatenate((beta, beta_p.reshape(*beta_p.shape, 1, 1)), axis=-2)

    # Return beta
    return(beta)

def get_epsilon(v,n,sigma):
    """
    Generates epsilon, the random error term.

    Parameters
    ----------
    v : int
        Number of voxels.
    n : int
        Number of observations.
    sigma : scalar sigma
        Standard deviation of the error term.

    Returns
    -------
    epsilon : numpy array for epsilon of dimensions (v,n,1).
    """

    # Get sigma2
    sigma2 = sigma**2

    # Make epsilon.
    epsilon = sigma2*np.random.randn(v,n)

    # Reshape to dimensions for broadcasting
    epsilon = epsilon.reshape(v, n, 1)

    return(epsilon)

def get_Y(X, beta, epsilon):
    """
    Generates the response vector Y.

    Parameters
    ----------
    X : numpy array for design matrix
    beta : numpy array for beta parameter vector
    epsilon : numpy array for epsilon of dimensions.

    Returns
    -------
    Y : numpy array for response vector
    """

    # Generate the response vector
    Y = X @ beta + epsilon

    # Return Y
    return(Y)

# Smoothing function
def smooth_data(data, D, fwhm, trunc=6, scaling='kernel'):

    # -----------------------------------------------------------------------
    # Reformat fwhm
    # -----------------------------------------------------------------------

    # Format fwhm and replace None with 0
    fwhm = np.asarray([fwhm]).ravel()
    fwhm = np.asarray([0. if elem is None else elem for elem in fwhm])

    # Non-zero dimensions
    D_nz = np.sum(fwhm>0)

    # Convert fwhm to sigma values
    sigma = fwhm / np.sqrt(8 * np.log(2))

    # -----------------------------------------------------------------------
    # Perform smoothing (this code is based on `_smooth_array` from the
    # nilearn package)
    # -----------------------------------------------------------------------
    
    # Loop through each dimension and smooth
    for n, s in enumerate(sigma):

        # If s is non-zero smooth by s in that direction.
        if s > 0.0:

            # Perform smoothing in nth dimension
            ndimage.gaussian_filter1d(data, s, output=data, mode='constant', axis=n, truncate=trunc)


    # -----------------------------------------------------------------------
    # Rescale
    # -----------------------------------------------------------------------
    if scaling=='kernel':
    
        # -----------------------------------------------------------------------
        # Rescale smoothed data to standard deviation 1 (this code is based on
        # `_gaussian_kernel1d` from the `scipy.ndimage` package).
        # -----------------------------------------------------------------------

        # Calculate sigma^2
        sigma2 = sigma*sigma

        # Calculate kernel radii
        radii = np.int16(trunc*sigma + 0.5)

        # Initialize array for phi values (has to be object as dimensions can 
        # vary in length)
        phis = np.empty(shape=(D_nz),dtype=object)

        # Index for non-zero dimensions
        j = 0

        # Loop through dimensions to get scaling constants
        for k in np.arange(D):

            # Skip the non-smoothed dimensions
            if fwhm[k]!=0:

                # Get range of values for this dimension
                r = np.arange(-radii[k], radii[k]+1)
                
                # Get the kernel for this dimension
                phi = np.exp(-0.5 / sigma2[k] * r ** 2)

                # Normalise phi
                phi = phi / phi.sum()

                # Add phi to dictionary
                phis[j]= phi[::-1]

                # Increment j
                j = j + 1
                
        # Create the D_nz dimensional grid
        grids = np.meshgrid(*phis);

        # Calculate the outer product of the kernels
        if D_nz == 1:

            # In case of only one dimension being smoothed, kernel_outer_product is simply the kernel itself.
            kernel_outer_product = phis[0]

            # Calculate the rescaling factor
            ss = np.sum(kernel_outer_product**2)

            # Rescale the data
            data = data / np.sqrt(ss)

        elif D_nz > 1:

            kernel_outer_product = np.outer(phis[0], phis[1])
            for j in range(2, D_nz):
                kernel_outer_product = np.outer(kernel_outer_product, phis[j])

            # Calculate the rescaling factor
            ss = np.sum(kernel_outer_product**2)

            # Rescale the data
            data = data / np.sqrt(ss)

    elif scaling=='max':

        # Rescale noise by dividing by maximum value
        data = data/np.max(data)

    return(data)