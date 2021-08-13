import numpy as np 
import os
import time
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from lib.generateData import *
from lib.fileio import *


# TODO:
# - cleanup

# field - field to thresh
# c - thresh
# d - dimension along which we get bdry
def get_bdry_map(field, c, d, mask=None): 

    # Get field dimensions
    dim = field.shape

    # Get boolean for where field is greater than c
    exc_set = field>c

    # ----------------------------------------------------------------------
    # Work out shift down and shift up indices
    # ----------------------------------------------------------------------
    # Initialize empty index collections
    up_indices = np.empty(shape=(len(dim)),dtype=object)
    down_indices = np.empty(shape=(len(dim)),dtype=object)


    # Loop through indices and join them in a list
    for i in np.arange(len(dim)):

        # If the d is not the i^th dimension keep all rows
        if d == i:

            # Add all but last to up indices
            up_indices[i] = np.arange(dim[i]-1)

            # Add all but first to down indices
            down_indices[i] = np.arange(1, dim[i])

        # If the d is not the i^th dimension keep all rows
        else:

            # Add full numpy range to both up and down
            up_indices[i] = np.arange(dim[i])
            down_indices[i] = np.arange(dim[i])

    # Padding dimensions to add zeros back on at the end (basically just 
    # same as dim but with the d^th dimension set to 1
    pdim = np.array(dim)
    pdim[d] = 1

    # ----------------------------------------------------------------------
    # Bottom boundary
    # ----------------------------------------------------------------------

    # Boundary without the bottom row (we can't assess downwards violations
    # when there is no data below)
    bdry_wo_bottom = exc_set[np.ix_(*up_indices)]

    # Complement of boundary downwards shifted
    shif_bdry_bottom = np.logical_not(exc_set[np.ix_(*down_indices)])

    # The boundary "below" the excursion set in this dimension
    bottom_bdry = bdry_wo_bottom*shif_bdry_bottom

    # If we have a mask we need to apply it
    if mask is not None:

        print('applying mask')
        mask = mask.reshape(field.shape)
        bottom_bdry = bottom_bdry*mask[np.ix_(*down_indices)]

    # Inner bottom bdry - add back on a bottom row
    bottom_bdry_inner = np.concatenate((bottom_bdry,np.zeros(pdim)),axis=d)

    # Outer bottom bdry - add back on a top row
    bottom_bdry_outer = np.concatenate((np.zeros(pdim),bottom_bdry),axis=d)

    # ----------------------------------------------------------------------
    # Top boundary
    # ----------------------------------------------------------------------

    # Boundary without the top row (we can't assess upwards violations
    # when there is no data above)
    bdry_wo_top = exc_set[np.ix_(*down_indices)]

    # Complement of boundary upwards shifted
    shif_bdry_top = np.logical_not(exc_set[np.ix_(*up_indices)])

    # The boundary "above" the excursion set in this dimension
    top_bdry = bdry_wo_top*shif_bdry_top

    # If we have a mask we need to apply it
    if mask is not None:
        mask = mask.reshape(field.shape)
        top_bdry = top_bdry*mask[np.ix_(*up_indices)]

    # Inner top bdry - add back on a top row
    top_bdry_inner = np.concatenate((np.zeros(pdim),top_bdry),axis=d)

    # Outer top bdry - add back on a bottom row
    top_bdry_outer = np.concatenate((top_bdry,np.zeros(pdim)),axis=d)

    # plt.figure(0)
    # plt.imshow(1*mask+2*top_bdry_inner+3*top_bdry_outer+4*bottom_bdry_inner+5*bottom_bdry_outer)
    # plt.savefig('/well/nichols/users/inf852/ConfSets2/tmp.png')
    # plt.figure(1)
    # plt.imshow(1*mask)
    # plt.savefig('/well/nichols/users/inf852/ConfSets2/tmp2.png')
    # plt.figure(2)
    # plt.imshow(2*top_bdry_inner+3*top_bdry_outer+4*bottom_bdry_inner+5*bottom_bdry_outer)
    # plt.savefig('/well/nichols/users/inf852/ConfSets2/tmp3.png')


    return(bottom_bdry_inner, bottom_bdry_outer, top_bdry_inner, top_bdry_outer)


def get_bdry_maps(field, c, mask=None):

    # Shape of field
    shape = np.array(field.shape)

    # Dimension of field
    dim = np.array(field.ndim)

    # Make an empty dictionary to store boundaries
    bdry_maps = dict()

    # Loop through dimensions of field and get the boundary boolean maps.
    for d in np.arange(dim):

        # Dimensions of 1 are assumed to be uninteresting as they are usually 
        # only included for broadcasting purposes.
        if shape[d]>1:

            # Get boundaries
            bottom_inner, bottom_outer, top_inner, top_outer = get_bdry_map(field, c, d, mask)

            # Record d^th boundary
            bdry_maps[d] = dict()

            # Add bottom boundaries
            bdry_maps[d]['bottom'] = dict()

            # Add inner and outer bottom boundaries
            bdry_maps[d]['bottom']['inner'] = bottom_inner
            bdry_maps[d]['bottom']['outer'] = bottom_outer

            # Add top boundaries
            bdry_maps[d]['top'] = dict()

            # Add inner and outer top boundaries
            bdry_maps[d]['top']['inner'] = top_inner
            bdry_maps[d]['top']['outer'] = top_outer

    # Add the non-flat (>1) dimensions as an array for good measure
    bdry_maps['dims'] = np.arange(dim)[shape>1]

    # Save original field shapes
    bdry_maps['shape_orig'] = np.array(field.shape)

    # Save original field dimension
    bdry_maps['dim_orig'] = np.array(field.ndim)

    # Return the bounaries
    return(bdry_maps)


def get_bdry_map_combined(field, c, mask=None):

    # Shape of field
    shape = np.array(field.shape)

    # Dimension of field
    dim = np.array(field.ndim)

    # Record if this is the first boundary we've seen
    first = True

    # Loop through dimensions of field and get the boundary boolean maps.
    for d in np.arange(dim):

        # Dimensions of 1 are assumed to be uninteresting as they are usually 
        # only included for broadcasting purposes.
        if shape[d]>1:

            # Get boundaries
            bottom_inner, bottom_outer, top_inner, top_outer = get_bdry_map(field, c, d, mask)

            # Record d^th boundary
            if first:

                # Add inner and outer bottom boundaries
                bdry_map = bottom_inner + bottom_outer + top_inner + top_outer

                # No longer looking at first boundary
                first = False

            else:

                # Add inner and outer bottom boundaries
                bdry_map = bdry_map + bottom_inner + bottom_outer + top_inner + top_outer

    # Return the boundary
    return(bdry_map>0)

def get_bdry_locs(bdry_maps):

    # Make an empty dictionary to store boundaries
    bdry_locs = dict()

    # Directions we can interpolate in
    directions = ['bottom', 'top']

    # -------------------------------------------------------------------------------------
    # Minimal dimensions (i.e. dimensions we would have if all dimensions of length 1 were
    # removed
    # -------------------------------------------------------------------------------------
    # Read in 0^th dimension to get a boundary image to work with
    dims = bdry_maps['dims']

    # Get the new shape we want to record indices for.
    new_shape = bdry_maps['shape_orig'][dims]

    # -------------------------------------------------------------------------------------
    # Loop through dimensions of field and get the locations of the boundary in the 
    # boundary boolean maps.
    # -------------------------------------------------------------------------------------
    for d in bdry_maps['dims']:

        # Record d^th boundary
        bdry_locs[d] = dict()

        # Loop through all directions getting locations
        for direction in directions:

            # Add bottom boundaries
            bdry_locs[d][direction] = dict()

            # Get inner map and reshape 
            inner = bdry_maps[d][direction]['inner'].reshape(new_shape)

            # Get outer map and reshape 
            outer = bdry_maps[d][direction]['outer'].reshape(new_shape)

            # Get coordinates of non-zero entries
            bdry_locs[d][direction]['inner'] = np.where(inner)
            bdry_locs[d][direction]['outer'] = np.where(outer)

    # Add the non-zero dimensions as an array for good measure
    bdry_locs['dims'] = bdry_maps['dims']

    # Return the bounaries
    return(bdry_locs)


def get_bdry_values(field, bdry_locs):

    # New dictionary to store boundary values
    bdry_vals = dict()

    # Directions we can interpolate in
    directions = ['bottom', 'top']

    # Boolean to tell if this is the first edge we are looking at.
    first = True

    # Loop through dimensions of field and get the boundary boolean maps.
    for d in bdry_locs['dims']:

        # Record d^th boundary
        bdry_vals[d] = dict()

        # Loop through all directions getting locations
        for direction in directions:

            # Record d^th boundary
            bdry_vals[d][direction] = dict()

            # Get inner and outer boundary values in this dimension and
            # direction            
            inner_vals = field[...,(*bdry_locs[d][direction]['inner'])]
            outer_vals = field[...,(*bdry_locs[d][direction]['outer'])]

            # Record boundary values
            bdry_vals[d][direction]['inner'] = inner_vals
            bdry_vals[d][direction]['outer'] = outer_vals

    # Add the non-zero dimensions as an array for good measure
    bdry_vals['dims'] = bdry_locs['dims']

    # Return boundary values
    return(bdry_vals)

def get_bdry_weights(bdry_vals,c):

    # New dictionary to store weights for interpolation along boundary
    bdry_weights = dict()

    # Directions we can interpolate in
    directions = ['bottom', 'top']

    # Loop through dimensions of field and get the boundary boolean maps.
    for d in bdry_vals['dims']:

        # Record d^th boundary
        bdry_weights[d] = dict()

        # Loop through all directions getting locations
        for direction in directions:

            # Record d^th boundary
            bdry_weights[d][direction] = dict()

            # Get inner and outer boundary values in this dimension and
            # direction            
            inner_vals = bdry_vals[d][direction]['inner']
            outer_vals = bdry_vals[d][direction]['outer']

            # Temporarily turn off divide by 0 warnings, we'll handle these below
            # (This is over cautious... it only really is a problem for plateua 
            # mu fields)
            with np.errstate(divide='ignore'):
                
                # Work out weights
                bdry_weights[d][direction]['outer']= (inner_vals-c)/(inner_vals-outer_vals)
                bdry_weights[d][direction]['inner']= (c-outer_vals)/(inner_vals-outer_vals)

            # In case we had 2 values which were the same (can happen when looking down
            # on ramp)
            inf_locs = np.isinf(bdry_weights[d][direction]['inner'])

            # Replace infs
            bdry_weights[d][direction]['inner'][inf_locs]=1
            bdry_weights[d][direction]['outer'][inf_locs]=0

    # Add the non-zero dimensions as an array for good measure
    bdry_weights['dims'] = bdry_vals['dims']

    # Return boundary values
    return(bdry_weights)



def get_bdry_vals_interpolated(bdry_vals,bdry_weights,dictform=False):

    # New dictionary to store weights for interpolation along boundary
    bdry_interp = dict()

    # Directions we can interpolate in
    directions = ['bottom', 'top']

    # Boolean to tell if this is the first edge we are looking at.
    first = True

    # Loop through dimensions of field and get the boundary boolean maps.
    for d in bdry_vals['dims']:

        # Record d^th boundary
        bdry_interp[d] = dict()

        # Loop through all directions getting locations
        for direction in directions:

            # Get inner and outer boundary values in this dimension and
            # direction            
            outer_vals = bdry_vals[d][direction]['outer']
            inner_vals = bdry_vals[d][direction]['inner']

            # Work out weights
            outer_weights = bdry_weights[d][direction]['outer']
            inner_weights = bdry_weights[d][direction]['inner']
        
            # Work out interpolated values
            bdry_interp[d][direction] = inner_weights*inner_vals + outer_weights*outer_vals

            # If this is the first edge we've looked at initialise 
            # concatenated interpolated boundary array
            if first:

                # Initialise array
                bdry_interp_concat = bdry_interp[d][direction]

                # We're no longer looking at the first edge
                first = False

            else:

                # Add the boundary values we just worked out
                bdry_interp_concat =  np.concatenate((bdry_interp_concat,bdry_interp[d][direction]),axis=-1)

    # If we have set dictform to true, return the boundry values in dictionary
    # form, i.e. preserving which edge each value came from. (more useful for
    # bug testing)
    if dictform:

        # Add the non-zero dimensions as an array for good measure
        bdry_interp['dims'] = bdry_vals['dims']

        # Return boundary values in dictionary form
        return(bdry_interp)

    # Otherwise just return them all in one big concatenated array (default -
    # more useful in practice)
    else:

        # Return boundary values in concatenated form
        return(bdry_interp_concat)



def get_bdry_values_concat(field, bdry_locs):

    # Directions we can interpolate in
    directions = ['bottom', 'top']

    # Boolean to tell if this is the first edge we are looking at.
    first = True

    # Loop through dimensions of field and get the boundary boolean maps.
    for d in bdry_locs['dims']:

        # Loop through all directions getting locations
        for direction in directions:

            # Get inner and outer boundary values in this dimension and
            # direction            
            inner_vals = field[...,(*bdry_locs[d][direction]['inner'])]
            outer_vals = field[...,(*bdry_locs[d][direction]['outer'])]

            # Reshape for concatenation
            inner_vals = inner_vals.reshape((*inner_vals.shape),1)
            outer_vals = outer_vals.reshape((*outer_vals.shape),1)
            
            # If this is the first edge we've looked at initialise 
            # concatenated interpolated boundary array
            if first:

                # Initialise array
                bdry_vals_concat = np.concatenate((inner_vals,outer_vals),axis=-1)

                # We're no longer looking at the first edge
                first = False

            else:

                # Initialise array
                current_bdry_vals = np.concatenate((inner_vals,outer_vals),axis=-1)

                # Add the boundary values we just worked out
                bdry_vals_concat =  np.concatenate((bdry_vals_concat,current_bdry_vals),axis=-2)

    # Return boundary values
    return(bdry_vals_concat)


def get_bdry_weights_concat(bdry_vals_concat,c):

    # Get inner and outer boundary values in this dimension and
    # direction            
    inner_vals = bdry_vals_concat[...,0]
    outer_vals = bdry_vals_concat[...,1]

    # Temporarily turn off divide by 0 warnings, we'll handle these below
    # (This is over cautious... it only really is a problem for plateua 
    # mu fields)
    with np.errstate(divide='ignore'):
        
        # Work out weights
        bdry_weights_inner_concat= (inner_vals-c)/(inner_vals-outer_vals)
        bdry_weights_outer_concat= (c-outer_vals)/(inner_vals-outer_vals)

        # Work out weights
        bdry_weights_inner_concat = bdry_weights_inner_concat.reshape((*bdry_weights_inner_concat.shape),1)
        bdry_weights_outer_concat = bdry_weights_outer_concat.reshape((*bdry_weights_outer_concat.shape),1)

    # In case we had 2 values which were the same (can happen when looking down
    # on ramp)
    inf_locs = np.isinf(bdry_weights_inner_concat)

    # Replace infs
    bdry_weights_inner_concat[inf_locs]=1
    bdry_weights_outer_concat[inf_locs]=0

    # Concatenate them
    bdry_weights_concat = np.concatenate((bdry_weights_outer_concat,bdry_weights_inner_concat),axis=-1) 

    # Return boundary values
    return(bdry_weights_concat)



def get_bdry_vals_interpolated_concat(bdry_vals_concat,bdry_weights_concat):

    # Get inner and outer boundary values in this dimension and
    # direction            
    outer_vals = bdry_vals_concat[...,0]
    inner_vals = bdry_vals_concat[...,1]

    # Work out weights
    outer_weights = bdry_weights_concat[...,0]
    inner_weights = bdry_weights_concat[...,1]
    

    # Work out interpolated values
    bdry_interp_concat = inner_weights*inner_vals + outer_weights*outer_vals

    # Return boundary values in concatenated form
    return(bdry_interp_concat)


def get_data_1field(muSpec,noiseSpec,dim):

    # Obtain the noise fields
    noise = get_noise(noiseSpec, dim)

    # Obtain mu
    mu = get_mu(muSpec, dim)
    
    # Create the data
    data = mu + noise

    # Return the data and mu
    return(data,mu)


















# Even circle points function, taken from:
# https://stackoverflow.com/questions/33510979/generator-of-evenly-spaced-points-in-a-circle-in-python
def circle_points(r, n):
    circles = []
    for r, n in zip(r, n):
        t = np.linspace(0, 2*np.pi, n, endpoint=False)
        x = np.round(r * np.cos(t))
        y = np.round(r * np.sin(t))
        circles.append(np.c_[x, y])
    return circles[0]








def testfn():
    
    # Set simulation number
    simNo = 17

    # Threshold
    c = 2

    # Output directory
    outdir = '/home/tommaullin/Documents/ConfRes'

    # Check if simulation directory exists, if not make it
    if not os.path.isdir(os.path.join(outdir, 'sim' + str(simNo))):
        os.mkdir(os.path.join(outdir, 'sim' + str(simNo)))

    # ----------------------------------------------------------------------------
    # Simulation 1 and 2
    # ----------------------------------------------------------------------------
    if simNo in [1,2]:

        # Set M
        m = 2

        # Check if figures directory exists, if not make it
        if not os.path.isdir(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample')):
            os.mkdir(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample'))

        # Create mu1 and mu2 specification
        mu1 = {}
        mu2 = {}

        # Create noise1 and noise2 specification
        noise1 = {}
        noise2 = {}

        # Add FWHM for noise
        noise1['FWHM'] = [0, 3, 3]

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = [0, 3, 3]

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Add mu1 type
        mu1['type'] = 'circle2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = np.array([5,5])

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center
        mu1['center'] = np.array([-15,0])

        # Add mu2 type
        mu2['type'] = 'circle2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = np.array([5,5])

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3

        # Add mu2 center 
        mu2['center'] = np.array([15,0])

        # Get fields 1 and 2
        data1, mu1 = get_data_1field(mu1, noise1, np.array([80,100,100]))
        data2, mu2 = get_data_1field(mu2, noise2, np.array([80,100,100]))

        # Get masks for fields 1 and 2
        mask1 = 1*(mu1[0,:,:]>c)
        mask2 = 1*(mu2[0,:,:]>c)

        # Save to files
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask1.csv'),mask1,remove=True)
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask2.csv'),mask2,remove=True)

        # Get masks for fields 1 and 2
        boundary1 = 1*get_bdry_map_combined(mu1, c)[0,:,:]
        boundary2 = 1*get_bdry_map_combined(mu2, c)[0,:,:]

        # Save to files
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry1.csv'),boundary1,remove=True)
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry2.csv'),boundary2,remove=True)

        # Get minimum field
        minfield = np.minimum(mu1,mu2)

        # Combined set
        maskAll = 1*(minfield[0,:,:]>c)

        # Save to file
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'maskall.csv'),maskAll,remove=True)

        # Combined boundary
        boundaryAll = 1*get_bdry_map_combined(minfield, c)[0,:,:]

        # Save to file
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdryall.csv'),boundaryAll,remove=True)

    # ----------------------------------------------------------------------------
    # Simulation 3 and 4
    # ----------------------------------------------------------------------------
    if simNo in [3,4]:

        # Set M
        m = 2

        # Check if figures directory exists, if not make it
        if not os.path.isdir(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample')):
            os.mkdir(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample'))

        # Create mu1 and mu2 specification
        mu1 = {}
        mu2 = {}

        # Create noise1 and noise2 specification
        noise1 = {}
        noise2 = {}

        # Add FWHM for noise
        noise1['FWHM'] = [0, 3, 3]

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = [0, 3, 3]

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = np.array([5,5])

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center (we only vary mu2)
        mu1['center'] = np.array([-15,0])

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = np.array([5,5])

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3

        # Add mu2 center 
        mu2['center'] = np.array([15,0])

        # Get fields 1 and 2
        data1, mu1 = get_data_1field(mu1, noise1, np.array([80,100,100]))
        data2, mu2 = get_data_1field(mu2, noise2, np.array([80,100,100]))

        # Get masks for fields 1 and 2
        mask1 = 1*(mu1[0,:,:]>c)
        mask2 = 1*(mu2[0,:,:]>c)

        # Save to files
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask1.csv'),mask1,remove=True)
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask2.csv'),mask2,remove=True)

        # Get masks for fields 1 and 2
        boundary1 = 1*get_bdry_map_combined(mu1, c)[0,:,:]
        boundary2 = 1*get_bdry_map_combined(mu2, c)[0,:,:]

        # Save to files
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry1.csv'),boundary1,remove=True)
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry2.csv'),boundary2,remove=True)

        # Get minimum field
        minfield = np.minimum(mu1,mu2)

        # Combined set
        maskAll = 1*(minfield[0,:,:]>c)

        # Save to file
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'maskall.csv'),maskAll,remove=True)

        # Combined boundary
        boundaryAll = 1*get_bdry_map_combined(minfield, c)[0,:,:]

        # Save to file
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdryall.csv'),boundaryAll,remove=True)

    # ----------------------------------------------------------------------------
    # Simulation 11 and 12
    # ----------------------------------------------------------------------------
    if simNo in [11,12]:

        # Set M
        m = 2

        # Check if figures directory exists, if not make it
        if not os.path.isdir(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample')):
            os.mkdir(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample'))

        # Create mu1 and mu2 specification
        mu1 = {}
        mu2 = {}

        # Create noise1 and noise2 specification
        noise1 = {}
        noise2 = {}

        # Gradient
        grad = 1

        # Add FWHM for noise
        noise1['FWHM'] = [0, 3, 3]

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = [0, 3, 3]

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Add mu1 type
        mu1['type'] = 'ramp2D'

        # Add mu1 fwhm
        mu1['orient'] = 'horizontal'

        # Add mu1 a
        mu1['a'] = '%.15f' % (2-grad)

        # Add mu1 b
        mu1['b'] = '%.15f' % (2+grad)

        # Add mu2 type
        mu2['type'] = 'ramp2D'

        # Add mu2 fwhm
        mu2['orient'] = 'vertical'

        # Add mu2 a
        mu2['a'] = '%.15f' % (2-grad)

        # Add mu2 b
        mu2['b'] = '%.15f' % (2+grad)

        # Get fields 1 and 2
        data1, mu1 = get_data_1field(mu1, noise1, np.array([80,100,100]))
        data2, mu2 = get_data_1field(mu2, noise2, np.array([80,100,100]))

        # Get masks for fields 1 and 2
        mask1 = 1*(mu1[0,:,:]>c)
        mask2 = 1*(mu2[0,:,:]>c)

        # Save to files
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask1.csv'),mask1,remove=True)
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask2.csv'),mask2,remove=True)

        # Get masks for fields 1 and 2
        boundary1 = 1*get_bdry_map_combined(mu1, c)[0,:,:]
        boundary2 = 1*get_bdry_map_combined(mu2, c)[0,:,:]

        # Save to files
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry1.csv'),boundary1,remove=True)
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry2.csv'),boundary2,remove=True)

        # Get minimum field
        minfield = np.minimum(mu1,mu2)

        # Combined set
        maskAll = 1*(minfield[0,:,:]>c)

        # Save to file
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'maskall.csv'),maskAll,remove=True)

        # Combined boundary
        boundaryAll = 1*get_bdry_map_combined(minfield, c)[0,:,:]

        # Save to file
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdryall.csv'),boundaryAll,remove=True)

    # ----------------------------------------------------------------------------
    # Simulation 15 and 16
    # ----------------------------------------------------------------------------
    if simNo in [15,16]:

        # Number of fields, m
        ms = np.arange(2,6)

        # Loop through all m values
        for m in ms:

            # ---------------------------------------------------------------
            # Mus
            # ---------------------------------------------------------------
            # Create empty specifications
            mus = {}

            # Loop through mus, adding each field in turn
            for i in np.arange(m):

                # Check if figures directory exists, if not make it
                if not os.path.isdir(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample')):
                    os.mkdir(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample'))

                # New empty dict
                mus['mu'+str(i+1)]={}

                # Mu type
                mus['mu'+str(i+1)]['type'] = 'circle2D' 

                # Mu FWHM
                mus['mu'+str(i+1)]['fwhm'] = np.array([5,5])

                # Mu r
                mus['mu'+str(i+1)]['r'] = 40

                # Mu magnitude
                mus['mu'+str(i+1)]['mag'] = 3

                # Get some evenly spaced center points
                centers = circle_points(np.array([25]),np.array([m]))

                # Mu center
                mus['mu'+str(i+1)]['center'] = centers[i,:]

            # ---------------------------------------------------------------
            # Epsilons
            # ---------------------------------------------------------------
            # Create empty specifications
            noises = {}

            # Loop through noises, adding each field in turn
            for i in np.arange(m):

                # New empty dict
                noises['noise'+str(i+1)]={}

                # Add FWHM
                noises['noise'+str(i+1)]['FWHM'] = [0, 3, 3]

                # Add type
                noises['noise'+str(i+1)]['type'] = 'homogen'

            # ---------------------------------------------------------------
            # Save images
            # ---------------------------------------------------------------
            if m == 2:

                # Get fields 1 and 2
                data1, mu1 = get_data_1field(mus['mu1'], noises['noise1'], np.array([80,100,100]))
                data2, mu2 = get_data_1field(mus['mu2'], noises['noise2'], np.array([80,100,100]))

                # Get masks for fields 1 and 2
                mask1 = 1*(mu1[0,:,:]>c)
                mask2 = 1*(mu2[0,:,:]>c)

                # Save to files
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask1.csv'),mask1,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask2.csv'),mask2,remove=True)

                # Get masks for fields 1 and 2
                boundary1 = 1*get_bdry_map_combined(mu1, c)[0,:,:]
                boundary2 = 1*get_bdry_map_combined(mu2, c)[0,:,:]

                # Save to files
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry1.csv'),boundary1,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry2.csv'),boundary2,remove=True)

                # Get minimum field
                minfield = np.minimum(mu1,mu2)

                # Combined set
                maskAll = 1*(minfield[0,:,:]>c)

                # Save to file
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'maskall.csv'),maskAll,remove=True)

                # Combined boundary
                boundaryAll = 1*get_bdry_map_combined(minfield, c)[0,:,:]

                # Save to file
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdryall.csv'),boundaryAll,remove=True)

            if m == 3:

                # Get fields 1 and 2
                data1, mu1 = get_data_1field(mus['mu1'], noises['noise1'], np.array([80,100,100]))
                data2, mu2 = get_data_1field(mus['mu2'], noises['noise2'], np.array([80,100,100]))
                data3, mu3 = get_data_1field(mus['mu3'], noises['noise3'], np.array([80,100,100]))

                # Get masks for fields 1 and 2
                mask1 = 1*(mu1[0,:,:]>c)
                mask2 = 1*(mu2[0,:,:]>c)
                mask3 = 1*(mu3[0,:,:]>c)

                # Save to files
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask1.csv'),mask1,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask2.csv'),mask2,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask3.csv'),mask3,remove=True)

                # Get masks for fields 1 and 2
                boundary1 = 1*get_bdry_map_combined(mu1, c)[0,:,:]
                boundary2 = 1*get_bdry_map_combined(mu2, c)[0,:,:]
                boundary3 = 1*get_bdry_map_combined(mu3, c)[0,:,:]

                # Save to files
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry1.csv'),boundary1,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry2.csv'),boundary2,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry3.csv'),boundary3,remove=True)

                # Get minimum field
                minfield = np.minimum(np.minimum(mu1,mu2),mu3)

                # Combined set
                maskAll = 1*(minfield[0,:,:]>c)

                # Save to file
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'maskall.csv'),maskAll,remove=True)

                # Combined boundary
                boundaryAll = 1*get_bdry_map_combined(minfield, c)[0,:,:]

                # Save to file
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdryall.csv'),boundaryAll,remove=True)

            if m == 4:

                # Get fields 1 and 2
                data1, mu1 = get_data_1field(mus['mu1'], noises['noise1'], np.array([80,100,100]))
                data2, mu2 = get_data_1field(mus['mu2'], noises['noise2'], np.array([80,100,100]))
                data3, mu3 = get_data_1field(mus['mu3'], noises['noise3'], np.array([80,100,100]))
                data4, mu4 = get_data_1field(mus['mu4'], noises['noise4'], np.array([80,100,100]))

                # Get masks for fields 1 and 2
                mask1 = 1*(mu1[0,:,:]>c)
                mask2 = 1*(mu2[0,:,:]>c)
                mask3 = 1*(mu3[0,:,:]>c)
                mask4 = 1*(mu4[0,:,:]>c)

                # Save to files
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask1.csv'),mask1,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask2.csv'),mask2,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask3.csv'),mask3,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask4.csv'),mask4,remove=True)

                # Get masks for fields 1 and 2
                boundary1 = 1*get_bdry_map_combined(mu1, c)[0,:,:]
                boundary2 = 1*get_bdry_map_combined(mu2, c)[0,:,:]
                boundary3 = 1*get_bdry_map_combined(mu3, c)[0,:,:]
                boundary4 = 1*get_bdry_map_combined(mu4, c)[0,:,:]

                # Save to files
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry1.csv'),boundary1,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry2.csv'),boundary2,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry3.csv'),boundary3,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry4.csv'),boundary4,remove=True)

                # Get minimum field
                minfield = np.minimum(np.minimum(np.minimum(mu1,mu2),mu3),mu4)

                # Combined set
                maskAll = 1*(minfield[0,:,:]>c)

                # Save to file
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'maskall.csv'),maskAll,remove=True)

                # Combined boundary
                boundaryAll = 1*get_bdry_map_combined(minfield, c)[0,:,:]

                # Save to file
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdryall.csv'),boundaryAll,remove=True)

            if m == 5:

                # Get fields 1 and 2
                data1, mu1 = get_data_1field(mus['mu1'], noises['noise1'], np.array([80,100,100]))
                data2, mu2 = get_data_1field(mus['mu2'], noises['noise2'], np.array([80,100,100]))
                data3, mu3 = get_data_1field(mus['mu3'], noises['noise3'], np.array([80,100,100]))
                data4, mu4 = get_data_1field(mus['mu4'], noises['noise4'], np.array([80,100,100]))
                data5, mu5 = get_data_1field(mus['mu5'], noises['noise5'], np.array([80,100,100]))

                # Get masks for fields 1 and 2
                mask1 = 1*(mu1[0,:,:]>c)
                mask2 = 1*(mu2[0,:,:]>c)
                mask3 = 1*(mu3[0,:,:]>c)
                mask4 = 1*(mu4[0,:,:]>c)
                mask5 = 1*(mu5[0,:,:]>c)

                # Save to files
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask1.csv'),mask1,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask2.csv'),mask2,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask3.csv'),mask3,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask4.csv'),mask4,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask5.csv'),mask5,remove=True)

                # Get masks for fields 1 and 2
                boundary1 = 1*get_bdry_map_combined(mu1, c)[0,:,:]
                boundary2 = 1*get_bdry_map_combined(mu2, c)[0,:,:]
                boundary3 = 1*get_bdry_map_combined(mu3, c)[0,:,:]
                boundary4 = 1*get_bdry_map_combined(mu4, c)[0,:,:]
                boundary5 = 1*get_bdry_map_combined(mu5, c)[0,:,:]

                # Save to files
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry1.csv'),boundary1,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry2.csv'),boundary2,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry3.csv'),boundary3,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry4.csv'),boundary4,remove=True)
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry5.csv'),boundary5,remove=True)

                # Get minimum field
                minfield = np.minimum(np.minimum(np.minimum(np.minimum(mu1,mu2),mu3),mu4),mu5)

                # Combined set
                maskAll = 1*(minfield[0,:,:]>c)

                # Save to file
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'maskall.csv'),maskAll,remove=True)

                # Combined boundary
                boundaryAll = 1*get_bdry_map_combined(minfield, c)[0,:,:]

                # Save to file
                append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdryall.csv'),boundaryAll,remove=True)

    # ----------------------------------------------------------------------------
    # Simulation 17 and 18
    # ----------------------------------------------------------------------------
    if simNo in [17,18]:

        # Set M
        m = 3

        # Check if figures directory exists, if not make it
        if not os.path.isdir(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample')):
            os.mkdir(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample'))

        # Create mu specification
        mus = {}

        # Create noise specification
        noises = {}

        # Empty mu spec
        mus['mu1'] = {}

        # Add mu1 type
        mus['mu1']['type'] = 'square2D' 

        # Add mu1 fwhm
        mus['mu1']['fwhm'] = np.array([5,5])

        # Add mu1 radius
        mus['mu1']['r'] = 30

        # Add mu1 magnitude
        mus['mu1']['mag'] = 3

        # Add mu1 center 
        mus['mu1']['center'] = np.array([-20,0])

        # Empty mu spec
        mus['mu2'] = {}

        # Add mu2 type
        mus['mu2']['type'] = 'square2D' 

        # Add mu2 fwhm
        mus['mu2']['fwhm'] = np.array([5,5])

        # Add mu2 radius
        mus['mu2']['r'] = 30

        # Add mu2 magnitude
        mus['mu2']['mag'] = 3

        # Add mu2 center 
        mus['mu2']['center'] = np.array([20,0])

        # Empty mu spec
        mus['mu3'] = {}

        # Add mu3 type
        mus['mu3']['type'] = 'square2D' 

        # Add mu3 fwhm
        mus['mu3']['fwhm'] = np.array([5,5])

        # Add mu3 radius
        mus['mu3']['r'] = 10

        # Add mu3 magnitude
        mus['mu3']['mag'] = 3

        # Add mu3 center 
        mus['mu3']['center'] = np.array([0,-20])

        # ---------------------------------------------------------------
        # Epsilons
        # ---------------------------------------------------------------
        # Create empty specifications
        noises = {}

        # Empty noise spec
        noises['noise1'] = {}

        # Add FWHM for noise
        noises['noise1']['FWHM'] = [0, 3, 3]

        # Add type for noise 1
        noises['noise1']['type'] = 'homogen'

        # Empty noise spec
        noises['noise2'] = {}

        # Add FWHM for noise
        noises['noise2']['FWHM'] = [0, 3, 3]

        # Add type for noise 1
        noises['noise2']['type'] = 'homogen'

        # Empty noise spec
        noises['noise3'] = {}

        # Add FWHM for noise
        noises['noise3']['FWHM'] = [0, 3, 3]

        # Add type for noise 1
        noises['noise3']['type'] = 'homogen'


        # Get fields 1 and 2
        data1, mu1 = get_data_1field(mus['mu1'], noises['noise1'], np.array([80,100,100]))
        data2, mu2 = get_data_1field(mus['mu2'], noises['noise2'], np.array([80,100,100]))
        data3, mu3 = get_data_1field(mus['mu3'], noises['noise3'], np.array([80,100,100]))

        # Get masks for fields 1 and 2
        mask1 = 1*(mu1[0,:,:]>c)
        mask2 = 1*(mu2[0,:,:]>c)
        mask3 = 1*(mu3[0,:,:]>c)

        # Save to files
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask1.csv'),mask1,remove=True)
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask2.csv'),mask2,remove=True)
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'mask3.csv'),mask3,remove=True)

        # Get masks for fields 1 and 2
        boundary1 = 1*get_bdry_map_combined(mu1, c)[0,:,:]
        boundary2 = 1*get_bdry_map_combined(mu2, c)[0,:,:]
        boundary3 = 1*get_bdry_map_combined(mu3, c)[0,:,:]

        # Save to files
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry1.csv'),boundary1,remove=True)
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry2.csv'),boundary2,remove=True)
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdry3.csv'),boundary3,remove=True)

        # Get minimum field
        minfield = np.minimum(np.minimum(mu1,mu2),mu3)

        # Combined set
        maskAll = 1*(minfield[0,:,:]>c)

        # Save to file
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'maskall.csv'),maskAll,remove=True)

        # Combined boundary
        boundaryAll = 1*get_bdry_map_combined(minfield, c)[0,:,:]

        # Save to file
        append_to_file(os.path.join(outdir, 'sim' + str(simNo), str(m) + 'sample', 'bdryall.csv'),boundaryAll,remove=True)

    # ----------------------------------------------------------------------------
    # Homogeneous noise
    # ----------------------------------------------------------------------------
    
    # Check if figures directory exists, if not make it
    if not os.path.isdir(os.path.join(outdir, 'homogen_noise')):
        os.mkdir(os.path.join(outdir, 'homogen_noise'))

    # Create noise specification
    noise = {}

    # Add FWHM for noise
    noise['FWHM'] = [0, 3, 3]

    # Add type for noise 
    noise['type'] = 'homogen'

    # Get some noise
    noise = get_noise(noise, np.array([80,100,100]))[0,:,:]

    # Save to files
    append_to_file(os.path.join(outdir, 'homogen_noise', 'data.csv'),noise,remove=True)

    # ----------------------------------------------------------------------------
    # Heterogeneous noise
    # ----------------------------------------------------------------------------
    
    # Check if figures directory exists, if not make it
    if not os.path.isdir(os.path.join(outdir, 'heterogen_noise')):
        os.mkdir(os.path.join(outdir, 'heterogen_noise'))

    # Create noise specification
    noise = {}

    # Add FWHM for noise
    noise['FWHM'] = [0, 3, 3]

    # Add type for noise 
    noise['type'] = 'heterogen'

    # Get some noise
    noise = get_noise(noise, np.array([80,100,100]))[0,:,:]

    # Save to files
    append_to_file(os.path.join(outdir, 'heterogen_noise', 'data.csv'),noise,remove=True)

    # -------------------------------------------------------------------
    # Color code:
    # -------------------------------------------------------------------
    # color = colornumber = combination of fields
    #
    # white = 0 = none
    # blue = 1 = 1
    # red = 2 = 2
    # orange = 3 = 3
    # green = 4 = 4
    # yellow = 5 = 5
    # purple = 6 = 1 + 2
    # sienna1 = 7 = 1 + 3
    # turquoise = 8 = 1 + 4
    # seagreen = 9 = 1 + 5
    # coral = 10 = 2 + 3
    # moccasin = 11 = 2 + 4
    # salmon = 12 = 2 + 5
    # rosybrown = 13 = 3 + 4
    # goldenrod = 14 = 3 + 5
    # khaki = 15 = 4 + 5
    # brown = 16 = all
    # -------------------------------------------------------------------


    # c = mcolors.ColorConverter().to_rgb
    # rvb = make_colormap([c('red'), c('orange'), c('yellow'), c('green'), c('blue'), c('violet'),c('red'), c('orange'), c('yellow'), c('green'), c('blue'), c('violet'),c('red'), c('orange'), c('yellow')])

    # plt.figure(0)
    # plt.imshow(demo_im,cmap=plt.get_cmap('tab20'))
    # plt.figure(1)
    # plt.imshow(demo_im,cmap=plt.get_cmap('tab20b'))
    # plt.figure(0)
    # plt.imshow(demo_im,cmap=plt.get_cmap('tab20c'))

    # plt.show()






    # mask = mu2[0,:,:]>c
 
    # bdryImage = get_bdry_map_combined(muHat1, c, mask)

    # plt.figure(0)
    # plt.imshow(bdryImage[0,:,:])

    # plt.figure(1)
    # plt.imshow(mask)

    # plt.figure(2)
    # plt.imshow(1*mask+1*bdryImage[0,:,:])

    # plt.show()

    # # -------------------------------------------------------------------

    # # Get the statistic field which defined Achat^{+/-,1/2}
    # stat1 = ((muHat1-c)/(sigma1*tau)).reshape(1,(*muHat1.shape))
    # stat2 = ((muHat2-c)/(sigma2*tau)).reshape(1,(*muHat2.shape))

    # # Minimum for intersection
    # stat = np.minimum(stat1,stat2) 
    # stat = stat.reshape(stat.shape[-2],stat.shape[-1])

    # # -------------------------------------------------------------------
    # # Get boolean maps for the boundary of Fc
    # Fc_bdry_maps = get_bdry_maps(np.minimum(mu1,mu2), c)

    # # Get coordinates for the boundary of Fc
    # Fc_bdry_locs = get_bdry_locs(Fc_bdry_maps)

    # # -------------------------------------------------------------------

    # # Obtain the values along the boundary for Fc
    # Fc_bdry_vals_concat = get_bdry_values_concat(np.minimum(mu1,mu2), Fc_bdry_locs)
    # #print('min mu vals concat: ', Fc_bdry_vals_concat)

    # # Obtain the weights along the boundary for Fc
    # Fc_bdry_weights_concat = get_bdry_weights_concat(Fc_bdry_vals_concat, c)
    # #print('min mu weights concat: ', Fc_bdry_weights_concat)

    # #print(stat.shape)

    # # Get the values along the outer and inner boundaries
    # stat_FcBdry1 = get_bdry_values_concat(stat, Fc_bdry_locs)
    # # #print('min stat vals concat: ', stat_FcBdry1)


    # # plt.figure(0)
    # # plt.hist(stat_FcBdry1.reshape(np.prod(stat_FcBdry1.shape)))#stat_FcBdry1.reshape(np.prod(stat_FcBdry.shape)))

    # #print('here2: ',stat_FcBdry1)

    # # Interpolate to get the values along the true boundary
    # stat_FcBdry1 = get_bdry_vals_interpolated_concat(stat_FcBdry1, Fc_bdry_weights_concat) # INCORRECT

    # # -------------------------------------------------------------------
    

    # # Obtain the values along the boundary for Fc
    # Fc_bdry_vals = get_bdry_values(np.minimum(mu1,mu2), Fc_bdry_locs)
    # # # print('min mu vals: ', Fc_bdry_vals)

    # # Obtain the weights along the boundary for Fc
    # Fc_bdry_weights = get_bdry_weights(Fc_bdry_vals, c)
    # # print('min mu weights: ', Fc_bdry_weights)
    
    # # Get the values along the outer and inner boundaries
    # stat_FcBdry1_tmp = get_bdry_values(stat, Fc_bdry_locs)
    # # print('min stat vals: ', stat_FcBdry1)






    # # Directions we can interpolate in
    # directions = ['bottom', 'top']

    # # Boolean to tell if this is the first edge we are looking at.
    # first = True

    # # Loop through dimensions of field and get the boundary boolean maps.
    # for d in Fc_bdry_locs['dims']:

    #     # Loop through all directions getting locations
    #     for direction in directions:

    #         stat_FcBdry1_tmp_concat_outer = stat_FcBdry1_tmp[d][direction]['outer']
    #         stat_FcBdry1_tmp_concat_inner = stat_FcBdry1_tmp[d][direction]['inner']

    #         stat_FcBdry1_tmp_concat_outer = stat_FcBdry1_tmp_concat_outer.reshape(1,stat_FcBdry1_tmp_concat_outer.shape[0])
    #         stat_FcBdry1_tmp_concat_inner = stat_FcBdry1_tmp_concat_inner.reshape(1,stat_FcBdry1_tmp_concat_inner.shape[0])

    #         if first:


    #             stat_FcBdry1_tmp_concat = np.concatenate((stat_FcBdry1_tmp_concat_outer,stat_FcBdry1_tmp_concat_inner),axis=0)

    #             first = False

    #         else:

    #             stat_FcBdry1_tmp_concat_current = np.concatenate((stat_FcBdry1_tmp_concat_outer,stat_FcBdry1_tmp_concat_inner),axis=0)
    #             stat_FcBdry1_tmp_concat = np.concatenate((stat_FcBdry1_tmp_concat_current,stat_FcBdry1_tmp_concat),axis=1)


    # # print('here: ',stat_FcBdry1_tmp_concat.transpose())

    # # plt.figure(0)
    # # plt.hist(stat_FcBdry1_tmp_concat.reshape(np.prod(stat_FcBdry1_tmp_concat.shape)))#stat_FcBdry2.reshape(np.prod(stat_FcBdry2.shape)))
    # # plt.show()



    # # Interpolate to get the values along the true boundary
    # stat_FcBdry1_tmp = get_bdry_vals_interpolated(stat_FcBdry1_tmp, Fc_bdry_weights) # CORRECT

    # # print(stat_FcBdry1.shape)
    # print(stat_FcBdry1_tmp.shape)

    # # print(stat_FcBdry1)
    # # print(stat_FcBdry1_tmp)

    # # print(np.all(stat_FcBdry1==stat_FcBdry1_tmp))

    # plt.figure(2)
    # plt.hist(stat_FcBdry1.reshape(np.prod(stat_FcBdry1.shape)))#stat_FcBdry2.reshape(np.prod(stat_FcBdry2.shape)))

    # plt.figure(3)
    # plt.hist(stat_FcBdry1_tmp.reshape(np.prod(stat_FcBdry1_tmp.shape)))#stat_FcBdry2.reshape(np.prod(stat_FcBdry2.shape)))
    # plt.show()
    # # # plt.figure(0)
    # # plt.imshow(muhat[0,:,:])

    # # plt.figure(1)
    # # plt.imshow(muhat[0,:,:]>2)

    # # plt.figure(2)
    # print('maps')
    # t1 = time.time()
    # bdry_maps = get_bdry_maps(muhat, c)
    # t2 = time.time()
    # print(t2-t1)

    # print('locs')
    # t1 = time.time()
    # bdry_locs = get_bdry_locs(bdry_maps)
    # t2 = time.time()
    # print(t2-t1)

    # print('vals')
    # t1 = time.time()
    # bdry_vals=get_bdry_values(muhat, bdry_locs)
    # t2 = time.time()
    # print(t2-t1)

    # print('weights')
    # t1 = time.time()
    # bdry_weights=get_bdry_weights(bdry_vals,c)
    # t2 = time.time()
    # print(t2-t1)

    # t1 = time.time()
    # bdry_vals=get_bdry_values(muhat, bdry_locs)
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

    # plt.show()

#testfn()