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
# - `fwhm`: Full Width Half Maximum for noise smoothness. Must be given as an
#           np array. Can include 0 or None for dimensions not to be
#           smoothed.
# - `dim`: Dimensions of data to be generated. Must be given as an np array.
#
# ===========================================================================
def get_data(muSpec,dim,fwhm):

    # Obtain the noise fields
    noise = get_noise(fwhm, dim)

    print(noise.shape)

    # Obtain mu
    mu = get_mu(muSpec, dim)
    
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
#              For 2D options, by default the returned field will be based 
#              on the last 2 dimensions of the `dim` argument.
# - `dim`: Dimensions of data to be generated. Must be given as an np array.
#
# ===========================================================================
def get_mu(muSpec, dim):

    # -----------------------------------------------------------------------
    # Generate mu
    # -----------------------------------------------------------------------

    # Ramp which increases from a to b uniformly
    if muSpec['type']=='ramp2D':

        # Get a and b
        a = muSpec['a']
        b = muSpec['b']

        # Ramp increasing vertically 
        if muSpec['orient']=='horizontal':

            # Using tile and transpose make the ramp
            mu = np.tile(np.linspace(a,b,dim[-2]).reshape(dim[-2],1),dim[-1]).transpose()

        # Ramp increasing horizontally
        elif muSpec['orient']=='vertical':

            # Using tile make the ramp
            mu = np.tile(np.linspace(a,b,dim[-1])[::-1].reshape(dim[-1],1),dim[-2])

    # Circular signal centered at 0 with radius 2.
    if muSpec['type']=='circle2D':

        # Radius
        r = muSpec['r']

        # Magnitude
        mag = muSpec['mag']

        # FWHM
        fwhm = muSpec['fwhm']

        # Adjust dimensions in case signal rolls off edge (we don't want signal to be smoothed
        # on the cutoff at the edge of the image as though there are zeros next to it)
        adjdim = np.array(dim)+2*r

        # Work out circle center (setting origin to the image center).
        center = np.array([adjdim[-2]//2, adjdim[-1]//2]) + muSpec['center']

        # Get an ogrid
        Y, X = np.ogrid[:adjdim[-2], :adjdim[-1]]

        # Make unsmoothed circular signal
        mu = np.array(np.sqrt((X-center[-2])**2+(Y-center[-1])**2) < r, dtype='float')

        # Smooth the data
        mu = mag*(smooth_data(mu, 2, fwhm, scaling='max'))

        # Recrop to original dimensions again
        mu = mu[...,r:adjdim[-2]-r,r:adjdim[-1]-r]


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
    trunc = 6

    # Format fwhm and replace None with 0
    fwhm = np.asarray([fwhm]).ravel()
    fwhm = np.asarray([0. if elem is None else elem for elem in fwhm])

    # -----------------------------------------------------------------------
    # Raw (padded) noise generation
    # -----------------------------------------------------------------------

    # Convert fwhm to sigma values
    sigma = fwhm / np.sqrt(8 * np.log(2))

    # Calculate kernel radii
    radii = np.int16(trunc*sigma + 0.5)

    # Work out padded dimensions
    pdim = dim + 2*(radii+1)

    # Generate unsmoothed random normal data for noise
    noise = np.random.randn(*pdim)

    # -----------------------------------------------------------------------
    # Perform smoothing
    # -----------------------------------------------------------------------
    noise = smooth_data(noise, D, fwhm, trunc)

    # -----------------------------------------------------------------------
    # Truncate the noise
    # -----------------------------------------------------------------------
    if D==2:
        noise = noise[(radii+1)[0]:(dim+radii+1)[0],(radii+1)[1]:(dim+radii+1)[1]]
    if D==3:
        noise = noise[(radii+1)[0]:(dim+radii+1)[0],(radii+1)[1]:(dim+radii+1)[1],(radii+1)[2]:(dim+radii+1)[2]]
    if D==4:
        noise = noise[(radii+1)[0]:(dim+radii+1)[0],(radii+1)[1]:(dim+radii+1)[1],(radii+1)[2]:(dim+radii+1)[2],(radii+1)[3]:(dim+radii+1)[3]]

    return(noise)

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

        # Initialize empty product grid
        product_grid = np.ones(grids[0].shape)

        # Loop through axes and take products
        for j in np.arange(D_nz):

            product_grid = grids[j]*product_grid

        # Get the normalizing constant by summing over grid
        ss = np.sum(product_grid**2)

        # Rescale noise
        data = data/np.sqrt(ss)

    elif scaling=='max':

        # Rescale noise by dividing by maximum value
        data = data/np.max(data)

    return(data)






def testfn():

    mu = get_mu({'type': 'circle2D_2mu', 'center': np.array([[-20,0],[20,0]]), 'fwhm': np.array([3,3]), 'r': np.array([25,25]), 'mag': np.array([3,3])}, np.array([1,100,100]))

    print(mu.shape)
    c=2

    plt.figure(0)
    plt.imshow(1*(mu[0,:,:]>c))
    plt.colorbar()
    
    plt.figure(1)
    plt.imshow(1*(mu[1,:,:]>c))
    plt.colorbar()
    
    plt.figure(2)
    plt.imshow(1*(mu[2,:,:]>c))
    plt.colorbar()
    
    plt.figure(3)
    plt.imshow(1*(mu[3,:,:]>c))
    plt.colorbar()

    plt.figure(4)
    plt.imshow(1*(mu[4,:,:]>c))
    plt.colorbar()

    plt.figure(5)
    plt.imshow(1*(mu[5,:,:]>c))
    plt.colorbar()

    # # plt.figure(2)
    # print('maps')
    # t1 = time.time()
    # bdry_maps = get_bdry_maps(mu, 2)
    # t2 = time.time()
    # print(t2-t1)

    # print('locs')
    # t1 = time.time()
    # bdry_locs = get_bdry_locs(bdry_maps)
    # t2 = time.time()
    # print(t2-t1)

    # print('vals')
    # t1 = time.time()
    # bdry_vals=get_bdry_values(mu, bdry_locs)
    # t2 = time.time()
    # print(t2-t1)

    # print('weights')
    # t1 = time.time()
    # bdry_weights=get_bdry_weights(bdry_vals,2)
    # t2 = time.time()
    # print(t2-t1)

    # t1 = time.time()
    # bdry_vals=get_bdry_values(mu, bdry_locs)
    # t2 = time.time()

    # print('interp')
    # t1 = time.time()
    # bdry_interp=get_bdry_vals_interpolated(bdry_vals,bdry_weights)
    # t2 = time.time()
    # print(t2-t1)

    # # print(bdry_maps)

    # # bdrys2 = get_bdry_maps(data, 2)

    # # print(bdrys2)
    # plt.figure(0)
    # plt.imshow(mu[0,:,:])
    # plt.colorbar()
    
    # plt.figure(1)
    # plt.imshow(muhat[0,:,:])
    # plt.colorbar()
    
    # plt.figure(2)
    # plt.imshow(data[8,:,:])
    # plt.colorbar()
    
    # plt.figure(3)
    # plt.imshow(mu[0,:,:]>2)
    
    # plt.figure(4)
    # plt.imshow(muhat[0,:,:]>2)

    # tmp = (bdry_maps[1]['top']['inner'][0,:,:]+bdry_maps[1]['bottom']['inner'][0,:,:]+bdry_maps[2]['top']['inner'][0,:,:]+bdry_maps[2]['bottom']['inner'][0,:,:])>0
    
    # plt.figure(5)
    # plt.imshow(tmp)

    # tmp2=(bdry_maps[1]['top']['outer'][0,:,:]+bdry_maps[1]['bottom']['outer'][0,:,:]+bdry_maps[2]['top']['outer'][0,:,:]+bdry_maps[2]['bottom']['outer'][0,:,:])>0
    # plt.figure(6)
    # plt.imshow(tmp2)

    # plt.figure(7)
    # plt.imshow(tmp+tmp2)

    plt.show()
