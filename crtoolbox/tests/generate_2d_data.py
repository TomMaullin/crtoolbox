import os
import numpy as np
from tests.generate_ni_data import smooth_data
from matplotlib import pyplot as plt

class Signal:
    """
    Base class for signal generation. This class is not meant to be used
    directly, but rather to be inherited by other classes.

    Attributes:
    -----------
    defaults : dict
        Dictionary containing default parameters for signal generation.
    params : dict
        Dictionary containing parameters for signal generation.
    mu : ndarray
        Array containing the generated signal.

    Methods:
    --------
    generate()
        Generate the signal.
    plot()
        Plot the generated signal.
    """
    
    def __init__(self, **kwargs):
        """
        Initialize the Signal class.
        
        Parameters:
        -----------
        kwargs : dict
            Dictionary containing parameters for signal generation.

        Returns:
        --------
        None
        """

        # Set the parameters
        self.params = self.defaults.copy()
        self.params.update(kwargs)

        # Initialize mu
        self.mu = None

        # Check if dim is 2D
        if len(self.params['dim']) > 2:
            raise ValueError("This class only supports 2D signal generation.")
    
    # Generate the signal
    def generate(self):
        """
        Generate the signal.

        Parameters:
        -----------
        None

        Returns:
        --------
        None
        """
        # Not implemented here
        raise NotImplementedError

    # Plot the signal
    def plot(self):
        """
        Plot the generated signal.

        Parameters:
        -----------
        None

        Returns:
        --------
        None
        """
        
        # Check if mu is generated
        if self.mu is None:
            # Generate mu
            self.generate()
        
        # Set the figure size
        plt.figure(figsize=(10, 6))

        # Plot the signal
        plt.imshow(self.mu)

        # Add a colorbar
        plt.colorbar()
        
        # Set the title and subtitle
        title = self.params['titlestr']
        subtitle = ', '.join([f"{k}={v}" for k, v in self.params.items() if k not in ['type', 'dim', 'titlestr']])
        
        # Add the title and subtitle
        plt.suptitle(title)
        plt.title(subtitle)

        # Show the plot
        plt.show()

# Circle signal class
class CircleSignal(Signal):
    """
    Class for generating a circle signal.
    
    Attributes:
    -----------
    defaults : dict
        Dictionary containing default parameters for signal generation.
    params : dict
        Dictionary containing parameters for signal generation.
    mu : ndarray
        Array containing the generated signal.
        
    Methods:
    --------
    generate()
        Generate the signal.
    plot()
        Plot the generated signal.

    Notes: 
    ------
    The following parameters are used for signal generation:
    - `r`: Radius of circle.
    - `mag`: Magnitude of signal.
    - `fwhm`: FWHM of smoothing kernel.
    - `center`: Center of circle.
    - `dim`: Dimensions of image.
    - `titlestr`: Title of plot.
    """

    # Default parameters
    defaults = {
        'type': 'circle2D', # Type of signal
        'r': 20, # Radius of circle
        'mag': 3, # Magnitude of signal
        'fwhm': np.array([6,6]), # fwhm of smoothing kernel
        'center': np.array([0, 0]), # Center of circle
        'dim': np.array([100, 100]), # Dimensions of image
        'titlestr': 'Circle Signal' # Title of plot
    }

    # Generate the signal
    def generate(self):
        """
        Generate the signal.

        Parameters:
        -----------
        None

        Returns:
        --------
        None

        Notes:
        ------
        The following parameters are used for signal generation:
        - `r`: Radius of circle.
        - `mag`: Magnitude of signal.
        - `fwhm`: FWHM of smoothing kernel.
        - `center`: Center of circle.
        - `dim`: Dimensions of image.
        """

        # Get the parameters
        r = self.params['r']
        mag = self.params['mag']
        fwhm = np.array(self.params['fwhm'])
        center = np.array(self.params['center'])
        dim = np.array(self.params['dim'])

        # Adjust dimensions in case signal rolls off edge (we don't want signal to be smoothed
        # on the cutoff at the edge of the image as though there are zeros next to it)
        adjdim = np.array(dim)+(np.array(dim)>1)*2*r

        # Work out circle center (setting origin to the image center).
        center = np.array([adjdim[-2]//2, adjdim[-1]//2]) + center

        # Get an ogrid
        Y, X = np.ogrid[:adjdim[-2], :adjdim[-1]]

        # Make unsmoothed circular signal
        mu = np.array(np.sqrt((X-center[-2])**2+(Y-center[-1])**2) < r, dtype='float')

        # Smooth the data
        mu = mag*(smooth_data(mu, 2, fwhm, scaling='max'))

        # Recrop to original dimensions again
        mu = mu[...,r:adjdim[-2]-r,r:adjdim[-1]-r]

        # Save mu
        self.mu = mu

# Square signal class
class SquareSignal(Signal):
    """
    Class for generating a square signal.

    Attributes:
    -----------
    defaults : dict
        Dictionary containing default parameters for signal generation.
    params : dict
        Dictionary containing parameters for signal generation.
    mu : ndarray
        Array containing the generated signal.

    Methods:
    --------
    generate()
        Generate the signal.
    plot()
        Plot the generated signal.
    
    Notes:
    ------
    The following parameters are used for signal generation:
    - `r`: Radius of square.
    - `mag`: Magnitude of signal.
    - `fwhm`: fwhm of smoothing kernel.
    - `center`: Center of square.
    - `dim`: Dimensions of image.
    - `titlestr`: Title of plot.
    """
    defaults = {
        'type': 'square2D',
        'r': 20,
        'mag': 3,
        'fwhm': np.array([6,6]),
        'center': np.array([0, 0]),
        'dim': np.array([100, 100]),
        'titlestr': 'Square Signal'
    }

    # Generate the signal
    def generate(self):
        """
        Generate the signal.

        Parameters:
        -----------
        None

        Returns:
        --------
        None

        Notes:
        ------
        The following parameters are used for signal generation:
        - `r`: Radius of square.
        - `mag`: Magnitude of signal.
        - `fwhm`: fwhm of smoothing kernel.
        - `center`: Center of square.
        - `dim`: Dimensions of image.
        """

        # Get the parameters
        r = self.params['r']
        mag = self.params['mag']
        fwhm = np.array(self.params['fwhm'])
        center = np.array(self.params['center'])
        dim = np.array(self.params['dim'])

        # Adjust dimensions in case signal rolls off edge (we don't want signal to be smoothed
        # on the cutoff at the edge of the image as though there are zeros next to it)
        adjdim = np.array(dim)+(np.array(dim)>1)*2*r

        # Work out square center (setting origin to the image center).
        center = np.array([adjdim[-2]//2, adjdim[-1]//2]) + center

        # Get an ogrid
        Y, X = np.ogrid[center[-2]-r:center[-2]+r, center[-1]-r:center[-1]+r]

        # Make unsmoothed square signal
        mu = np.zeros((adjdim[-2],adjdim[-1]))
        mu[...,X,Y] = 1

        # Smooth the data and multiply by magnitude
        mu = mag*(smooth_data(mu, 2, fwhm, scaling='max'))

        # Recrop to original dimensions again
        mu = mu[...,r:adjdim[-2]-r,r:adjdim[-1]-r]

        # Save mu
        self.mu = mu

# Ramp signal class
class RampSignal(Signal):
    """
    Class for generating a ramp signal.

    Attributes:
    -----------
    defaults : dict
        Dictionary containing default parameters for signal generation.
    params : dict
        Dictionary containing parameters for signal generation.
    mu : ndarray
        Array containing the generated signal.

    Methods:
    --------
    generate()
        Generate the signal.
    plot()
        Plot the generated signal.

    Notes:
    ------
    The following parameters are used for signal generation:
    - `a`: Lower value of ramp.
    - `b`: Upper value of ramp.
    - `orient`: Orientation of ramp.
    - `dim`: Dimensions of image.
    - `titlestr`: Title of plot.
    """

    # Default parameters
    defaults = {
        'type': 'ramp2D', # Type of signal
        'a': 1, # Lower value of ramp
        'b': 3, # Upper value of ramp
        'orient': 'horizontal', # Orientation of ramp
        'dim': np.array([100, 100]), # Dimensions of image
        'titlestr': 'Ramp Signal' # Title of plot
    }

    # Generate the signal
    def generate(self):

        # Get the parameters
        a = float(self.params['a'])
        b = float(self.params['b'])
        dim = self.params['dim']
        
        # Ramp increasing vertically
        if self.params['orient'] == 'horizontal':
            mu = np.tile(np.linspace(a,b,dim[-2]).reshape(dim[-2],1),dim[-1]).transpose()
        # Ramp increasing horizontally
        elif self.params['orient'] == 'vertical':
            mu = np.tile(np.linspace(a,b,dim[-1])[::-1].reshape(dim[-1],1),dim[-2])

        # Save mu
        self.mu = mu

# Noise class
class Noise:
    """
    Class for noise generation.

    Attributes:
    -----------
    defaults : dict
        Dictionary containing default parameters for noise generation.
    params : dict
        Dictionary containing parameters for noise generation.
    noise : ndarray
        Array containing either the generated noise, or one instance of the
        generated noise.

    Methods:
    --------
    generate()
        Generate the noise, or one instance of the noise.
    plot()
        Plot an instance of the generated noise.

    Notes:
    ------
    The following parameters are used for noise generation:
    - `type`: Type of noise to generate. Can be `homogen` or `heterogen`.
    - `fwhm`: fwhm of smoothing kernel.
    - `var`: Variance of noise.
    - `dim`: Dimensions of image.
    - `titlestr`: Title of plot.
    """

    # Default parameters
    defaults = {
        'type': 'homogen',
        'fwhm': np.array([5,5]),
        'var': 5,
        'dim': np.array([100, 100]),
        'titlestr': 'Homogeneous Noise',
        'n': 100
    }

    # Initialize the noise
    def __init__(self, **kwargs):
        """
        Initialize the Noise class.

        Parameters:
        -----------
        kwargs : dict
            Dictionary containing parameters for noise generation.

        Returns:
        --------
        None
        """

        # Set the parameters
        self.params = self.defaults.copy()
        self.params.update(kwargs)

        # Initialize noise
        self.noise = None

    # Generate the noise
    def generate(self, single_slice=False):
        """
        Generate the noise, or one instance of the noise.

        Parameters:
        -----------
        single_slice : bool
            Whether to generate one instance of the noise, or the noise
            for all n observations.

        Returns:
        --------
        None

        Notes:
        ------
        The following parameters are used for noise generation:
        - `type`: Type of noise to generate. Can be `homogen` or `heterogen`.
        - `fwhm`: fwhm of smoothing kernel.
        - `var`: Variance of noise.
        - `dim`: Dimensions of image.
        - `titlestr`: Title of plot.
        """

        # Get the parameters
        fwhm = np.array(self.params['fwhm'])
        var = np.array(self.params['var'])
        dim = np.array(self.params['dim'])

        # If single_slice is True then generate one instance of the noise
        if single_slice:

            # Set n to 1
            n = np.array([1])

        # Otherwise generate n instances of the noise
        else:

            # Get n
            n = np.array(self.params['n']).reshape(1)

        # Concenate n and dim
        dim = np.concatenate((n, dim))

        # Record zero fwhm for dimensions not to be smoothed
        fwhm = np.concatenate((np.array([0]), fwhm))

        # Get dimension of data
        D = np.prod(np.shape(dim))

        # Truncation (how many standard deviations will the Gaussian filter be
        # truncated at)
        trunc = 6

        # Format fwhm and replace None with 0
        fwhm = np.asarray([fwhm]).ravel()
        fwhm = np.asarray([0. if elem is None else elem for elem in fwhm])

        # Convert fwhm to sigma values
        sigma = fwhm / np.sqrt(8 * np.log(2))

        # Calculate kernel radii
        radii = np.array(2*np.round_(trunc*sigma) + 1, dtype=int)

        # Work out padded dimensions
        pdim = dim + 2*(radii+1)

        # Generate unsmoothed random normal data for noise
        noise = np.random.randn(*pdim)

        # Perform smoothing
        noise = smooth_data(noise, D, fwhm, trunc)

        # Truncate the noise
        slice_indices = [slice((radii[d] + 1), (dim[d] + radii[d] + 1)) for d in range(D)]
        noise = noise[tuple(slice_indices)]

        # Heterogenous noise (ramp)
        if self.params['type'] == 'heterogen':

            # Multiply by a ramp
            noise = noise*np.linspace(0.5,1.5,noise.shape[-1])

            # Change title
            self.params['titlestr'] = 'Heterogeneous Noise'

        # Alter the variance of the noise
        noise = noise*var

        # Save the noise
        self.noise = noise

    # Plot the noise
    def plot(self):
        """
        Plot an instance of the generated noise.

        Parameters:
        -----------
        None

        Returns:
        --------
        None

        Notes:
        ------
        The following parameters are used for noise generation:
        - `type`: Type of noise to generate. Can be `homogen` or `heterogen`.
        - `fwhm`: fwhm of smoothing kernel.
        - `var`: Variance of noise.
        - `dim`: Dimensions of image.
        - `titlestr`: Title of plot.
        """

        # Check if noise is generated
        if self.noise is None:

            # Generate a noise instance if it isn't
            self.generate(single_slice=True)

        # Set the figure size
        plt.figure(figsize=(10, 6))

        # Plot the first slice
        plt.imshow(self.noise[0,...])  

        # Add a colorbar
        plt.colorbar()

        # Set the title and subtitle
        title = self.params['titlestr']
        subtitle = ', '.join([f"{k}={v}" for k, v in self.params.items() if k not in ['titlestr', 'dim']])
        
        # Add the title and subtitle
        plt.suptitle(title)
        plt.title(subtitle)

        # Show the plot
        plt.show()


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


# Generate 2D data
def generate_data_2D(mu, noise, out_dir=None):
    """
    Generate 2D data.

    Inputs:
    -------
    - `mu`: class object specifying the mean field.
    - `noise`: class object specifying the noise field.

    Outputs:
    --------
    - `data_files`: List of filenames for the generated data.
    - `mu_file`: Filename for the generated mean field.
    """

    # If out_dir is not given then set it to the current directory
    if out_dir is None:
        out_dir = os.getcwd()
    
    # Check if out_dir exists
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Obtain the noise fields
    if noise.noise is not None:

        # Check the first dimension of noise
        if noise.noise.shape[0] == 1:

            # Generate the noise
            noise.generate()

        # Get the noise
        noise = noise.noise

    else:
    
        # Generate the noise
        noise.generate()

        # Get the noise
        noise = noise.noise

    # Obtain mu
    if mu.mu is None:

        # Generate mu
        mu.generate()

    # Get mu
    mu = mu.mu

    # Append an extra dimension on the front of mu
    mu = mu.reshape((1,) + mu.shape)

    # Create the data
    data = mu + noise

    # Empty list for filenames
    data_files = []

    # Loop through first dimension and save to npy files
    for i in np.arange(data.shape[0]):

        # Save the data
        np.save(os.path.join(out_dir,'data'+str(i)+'.npy'),data[i,:,:])

        # Append the filename to the list
        data_files.append(os.path.join(out_dir,'data'+str(i)+'.npy'))

    # Save the mu
    np.save(os.path.join(out_dir,'mu.npy'),mu)

    # Record the mu filename
    mu_file = os.path.join(out_dir,'mu.npy')

    # Return the data and mu
    return(data_files,mu_file)

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
def get_data_2field(muSpec1,muSpec2,noiseSpec1,noiseSpec2,dim,noiseCorr=None):

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

    # Return the noises
    return(new_noise1,new_noise2)
