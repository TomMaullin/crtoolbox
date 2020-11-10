import numpy as np
from scipy import ndimage
import time
from matplotlib import pyplot as plt

# ===========================================================================
#
# Inputs:
#
# ---------------------------------------------------------------------------
#
# - `fwhm`: Full Width Half Maximum for noise smoothness. Must be given as an
#           np array. Can include 0 or None for dimensions not to be
#           smoothed.
# - `dim`: Dimensions of data to be generated. Must be given as an np array.
#
# ===========================================================================
def get_data(muType,dim,fwhm):

    # Obtain the noise fields
    noise = get_noise(fwhm, dim)

    # Obtain mu
    mu = get_mu(muType, dim)
    
    # Create the data
    data = mu + noise

    # Return the data and mu
    return(data, mu)


# ===========================================================================
#
# Inputs:
#
# ---------------------------------------------------------------------------
#
# - `muType`: Type of signal to generate. Current supported options are:
#              - `rampHoriz2D`: A linear ramp, increasing across the first
#                               dimension of the 2D mu field.
#              - `rampVert2D`: A linear ramp, increasing across the second
#                              dimension of the mu field.
#              For 2D options, by default the returned field will be based 
#              on the last 2 dimensions of the `dim` argument.
# - `dim`: Dimensions of data to be generated. Must be given as an np array.
#
# ===========================================================================
def get_mu(muType, dim):

    # -----------------------------------------------------------------------
    # Generate mu
    # -----------------------------------------------------------------------

    # Ramp increasing vertically from 1 to 3 horizontally
    if muType=='rampHoriz2D':

        # Using tile and transpose make the ramp
        mu = np.tile(np.linspace(1,3,dim[-2]).reshape(dim[-2],1),dim[-1]).transpose()

    # Ramp increasing vertically from 1 to 3 vertically
    elif muType=='rampVert2D':

        mu = np.tile(np.linspace(1,3,dim[-1])[::-1].reshape(dim[-1],1),dim[-2])

    # -----------------------------------------------------------------------
    # Give mu the correct dimensions to be broadcasted with the data we are
    # creating
    # -----------------------------------------------------------------------
    # Get the dimension of the mu field
    muD = mu.ndim

    # Get the dimension of the data we're creating
    dataD = np.prod(dim.shape)

    # Check what shape a needs to be to be output.
    if muD < dataD:

        # Work out new dimensions for image
        newD = np.ones(dataD,dtype=np.int)

        # Replace last few dimensions with mu dimensions
        for i in np.arange(muD):

            # Replace ith dimension
            newD[-(i+1)] = mu.shape[-(i+1)]

    # Reshape mu 
    mu = mu.reshape(newD)

    # Return mu
    return(mu)

# ===========================================================================
#
# Inputs:
#
# ---------------------------------------------------------------------------
#
# - `fwhm`: Full Width Half Maximum for noise smoothness. Must be given as an
#           np array. Can include 0 or None for dimensions not to be
#           smoothed.
# - `dim`: Dimensions of data to be generated. Must be given as an np array.
#
# ===========================================================================
def get_noise(fwhm, dim):

    # -----------------------------------------------------------------------
    # Useful scalars
    # -----------------------------------------------------------------------
    # Dimension D
    D = np.prod(np.shape(dim))

    # Truncation (how many standard deviations will the Gaussian filter be 
    # truncated at)
    trunc = 4

    # Format fwhm and replace None with 0
    fwhm = np.asarray([fwhm]).ravel()
    fwhm = np.asarray([0. if elem is None else elem for elem in fwhm])

    # Non-zero dimensions
    D_nz = np.sum(fwhm>0)

    # Convert fwhm to sigma values
    sigma = fwhm / np.sqrt(8 * np.log(2))

    # Calculate kernel radii
    radii = np.int16(trunc*sigma + 0.5)
    
    # -----------------------------------------------------------------------
    # Raw noise generation
    # -----------------------------------------------------------------------
    
    # Work out padded dimensions
    pdim = dim + 2*(radii+1)

    # Generate unsmoothed random normal data for noise
    noise = np.random.randn(*pdim)

    # -----------------------------------------------------------------------
    # Perform smoothing (this code is based on `_smooth_array` from the
    # nilearn package)
    # -----------------------------------------------------------------------
    
    # Loop through each dimension and smooth
    for n, s in enumerate(sigma):

        # If s is non-zero smooth by s in that direction.
        if s > 0.0:

            # Perform smoothing in nth dimension
            ndimage.gaussian_filter1d(noise, s, output=noise, mode='constant', axis=n, truncate=trunc)

    # Truncate the noise
    if D==2:
        noise = noise[(radii+1)[0]:(dim+radii+1)[0],(radii+1)[1]:(dim+radii+1)[1]]
    if D==3:
        noise = noise[(radii+1)[0]:(dim+radii+1)[0],(radii+1)[1]:(dim+radii+1)[1],(radii+1)[2]:(dim+radii+1)[2]]
    if D==4:
        noise = noise[(radii+1)[0]:(dim+radii+1)[0],(radii+1)[1]:(dim+radii+1)[1],(radii+1)[2]:(dim+radii+1)[2],(radii+1)[3]:(dim+radii+1)[3]]

    # -----------------------------------------------------------------------
    # Rescale smoothed data to standard deviation 1 (this code is based on
    # `_gaussian_kernel1d` from the `scipy.ndimage` package.
    # -----------------------------------------------------------------------

    # Calculate sigma^2
    sigma2 = sigma*sigma

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

    # Initialize empty product grid
    product_grid = np.ones(grids[0].shape)

    # Loop through axes and take products
    for j in np.arange(D_nz):

        product_grid = grids[j]*product_grid

    # Get the normalizing constant by summing over grid
    ss = np.sum(product_grid**2);

    # Rescale noise
    noise = noise/np.sqrt(ss)

    return(noise)
