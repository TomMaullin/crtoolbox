import numpy as np
from scipy import ndimage
import time
from matplotlib import pyplot as plt
from tests.generate_ni_data import smooth_data



"""
Even circle points function, taken from:
https://stackoverflow.com/questions/33510979/generator-of-evenly-spaced-points-in-a-circle-in-python

Inputs:
- `r`: Radius of circle.
- `n`: Number of points to generate.

Outputs:
- `circles`: Array of circle points.
"""
def circle_points(r, n):
    circles = []
    for r, n in zip(r, n):
        t = np.linspace(0, 2*np.pi, n, endpoint=False)
        x = np.round(r * np.cos(t))
        y = np.round(r * np.sin(t))
        circles.append(np.c_[x, y])
    return circles[0]


def get_data_1field(muSpec,noiseSpec,dim):

    # Obtain the noise fields
    noise = get_noise(noiseSpec, dim)

    # Obtain mu
    mu = get_mu(muSpec, dim)
    
    # Create the data
    data = mu + noise

    # Return the data and mu
    return(data,mu)

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
def get_data(muSpec1,muSpec2,noiseSpec1,noiseSpec2,dim,noiseCorr=None):

    # Obtain the noise fields
    noise1 = get_noise(noiseSpec1, dim)
    noise2 = get_noise(noiseSpec2, dim)

    # Correlate the data if needed
    if noiseCorr is not None:
        noise1, noise2 = correlateData(noise1,noise2,noiseCorr)

    # Obtain mu
    mu1 = get_mu(muSpec1, dim)
    mu2 = get_mu(muSpec2, dim)
    
    # Create the data
    data1 = mu1 + noise1
    data2 = mu2 + noise2

    # Return the data and mu
    return(data1, data2, mu1, mu2)


"""
Function to get mean field mu from a specification dictionary.

Inputs:
- `muSpec`: Dictionary specifying mu to simulate. Always must include a
            `type` parameter specifying `ramp2D` or `circle2D`.
- `dim`: Dimensions of data to be generated. Must be given as an np array.

Outputs:
- `mu`: Mean field mu.
"""
def get_mu(muSpec, dim):

    # -----------------------------------------------------------------------
    # Generate mu
    # -----------------------------------------------------------------------

    # Ramp which increases from a to b uniformly
    if muSpec['type']=='ramp2D':

        # Get a and b
        a = np.float(muSpec['a'])
        b = np.float(muSpec['b'])

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

    # Square signal centered at 0 with radius 2.
    if muSpec['type']=='square2D':

        # Radius
        r = muSpec['r']

        # Magnitude
        mag = muSpec['mag']

        # FWHM
        fwhm = muSpec['fwhm']

        # Adjust dimensions in case signal rolls off edge (we don't want signal to be smoothed
        # on the cutoff at the edge of the image as though there are zeros next to it)
        adjdim = np.array(dim)+(np.array(dim)>1)*2*r

        # Work out square center (setting origin to the image center).
        center = np.array([adjdim[-2]//2, adjdim[-1]//2]) + muSpec['center']

        # Get an ogrid
        Y, X = np.ogrid[center[-2]-r:center[-2]+r, center[-1]-r:center[-1]+r]

        # Make unsmoothed square signal
        mu = np.zeros((adjdim[-2],adjdim[-1]))
        mu[...,X,Y]=1

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

"""
Function to generate noise based on a given noise specification.

Inputs:
- `noiseSpec`: Dictionary specifying noise to simulate. Always must include:
    -  a `type` parameter specifying `heterogen` or `homogen` for 
       heterogenous or homogenous noise respectively. 
    -  a `FWHM` parameter specifying the FWHM of the Gaussian filter used to 
       smooth the noise. If `None` is given then no smoothing is applied.
    -  a `mag` parameter specifying the magnitude of the noise.
- `dim`: Dimensions of data to be generated. Must be given as an np array.

Outputs:
- `noise`: Noise field generated based on the given noise specification.
"""
def get_noise(noiseSpec, dim):

    # Get FWHM
    fwhm = noiseSpec['FWHM']

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
    radii = np.array(2*np.round_(trunc*sigma) + 1,dtype=int)

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

    # Heterogenous noise (ramp)
    if noiseSpec['type']=='heterogen':

        noise = noise*np.linspace(0.5,1.5,noise.shape[-1])

    # Alter the magnitude of the noise
    if 'mag' in noiseSpec:
        noise = noise*np.float(noiseSpec['mag'])

    return(noise)

"""
Function to induce correlation between two noise fields.

Inputs:
- `noise1`: First noise field to correlate.
- `noise2`: Second noise field to correlate.
- `noiseCorr`: Correlation between the two noise fields.

Outputs:
- `new_noise1`: First noise field after correlation.
- `new_noise2`: Second noise field after correlation.
"""
def correlateData(noise1,noise2,noiseCorr):

    # Reshape noises
    noise1 = noise1.reshape(*noise1.shape,1)
    noise2 = noise2.reshape(*noise2.shape,1)

    # Combine them (for broadcasting covariance matrix multiplication)
    noises = np.concatenate((noise1,noise2),axis=-1)
    noises = noises.reshape(*noises.shape,1)

    abs_noiseCorr = np.abs(noiseCorr)
    sgn_noiseCorr = np.sign(noiseCorr)

    # Work out covariance matrix we need.
    covMat = np.array([[1,0],[abs_noiseCorr, np.sqrt(1/abs_noiseCorr**2 -1)*abs_noiseCorr]])

    # Correlate noises
    new_noises = covMat @ noises

    # Get back noise 1 and noise 2
    new_noise1 = new_noises[...,0,0]
    new_noise2 = new_noises[...,1,0]

    # Multiply one of the fields by -1 if need negative correlation
    if (sgn_noiseCorr == -1):
        new_noise2 = -new_noise2
        print('-ve')

    # Return the noises
    return(new_noise1,new_noise2)
