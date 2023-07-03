import numpy as np 
import os
import time
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from crtoolbox.generate import *
from crtoolbox.lib.fileio import *

def get_bdry_map(field, c, d, mask=None): 
    """
    Given a field, threshold and dimension, the below function derives pairs of
    neighbouring pixels which lie inside and outside of the excursion set given
    by {pixels in field: field at pixel > c}.

    Parameters: 
    -----------
    - field: An image of a field  from which we wish to obtain an image of the
            boundary of the excursion set (for example, if we want the boundary 
            set {pixels : mean at pixel = c}, then `field' should be an image
            of the mean.
    - c: The value c with which to threshold the field (e.g. we are interested
        in the set {pixels : field at pixel = c}.
    - d: The dimension along which to derive boundary values (for example, in 
        2D if dimension = 0, this function will return horizontal pairs of
        neighbouring pixels which lie inside and outside the boundary
        respectively, whereas if dimension = 1, the returned pairs will be
        vertical). Note: d is zero indexed.
    - mask: A binary image showing pixels which should be ignored during
            computation, or ``masked out". For example, if we are interested in
            fMRI BOLD signal over the brain, the mask might be only the pixels
            which lie inside the brain, as our signal is meaningless for pixels
            outside this region.

    Returns:
    --------
    - Four binary images 1-0 images of;
    - bottom_bdry_inner: The pixels which lie inside the excursion set at the
                            `bottom' in dimension d (e.g. the points labelled `3'
                            in the below image).
    - bottom_bdry_outer: The pixels which lie outside the excursion set at the
                            `bottom' in dimension d (e.g. the points labelled `4'
                            in the below image).
    - top_bdry_inner: The pixels which lie inside the excursion set at the
                        `top' in dimension d (e.g. the points labelled `2' in 
                        the below image).
    - top_bdry_outer: The pixels which lie outside the excursion set at the
                        `top' in dimension d (e.g. the points labelled `1' in 
                        the below image).
        _______________________________________________________________
        |                                                               |    
        |               11111111111111111111111                         |   | 
        |            111222222222222222222222221                        |   |  
        |           1222                       21                       |   |
        |           2                           2111                    |   |
        |          1                             2221                   |   |
        |          2                                21                  |   |
        |         |                                  21                 |   |
        |         \                                   2|                |   |  d
        |          |               A_c                 |                |   |
        |          |3                                  |                |   |
        |           43                                 3                |   |
        |            43                               34                |   |
        |             43                             34                 |   |
        |              43                           34                  |   |
        |               433                      3334                   |   |
        |                443                    3444                    |  \|/
        |                  4333333333333333333334                       |   V
        |                   44444444444444444444                        |
        |_______________________________________________________________|   
        Fig A. An image of the boundary values returned along dimension d for
        an excursion set, A_c.
    """

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

    return(bottom_bdry_inner, bottom_bdry_outer, top_bdry_inner, top_bdry_outer)

def get_bdry_maps(field, c, mask=None):
    """
    Given a field and a threshold, the below function returns all pairs of 
    neighbouring pixels which lie inside and outside of the excursion set given
    by {pixels in field: field at pixel > c}, respectively.

    Parameters:
    ----------
    - field: An image of a field  from which we wish to obtain an image of the
            boundary of the excursion set (for example, if we want the boundary 
            set {pixels : mean at pixel = c}, then `field' should be an image
            of the mean.
    - c: The value c with which to threshold the field (e.g. we are interested
        in the set {pixels : field at pixel = c}.

    - mask: A binary image showing pixels which should be ignored during
            computation, or ``masked out". For example, if we are interested in
            fMRI BOLD signal over the brain, the mask might be only the pixels
            which lie inside the brain, as our signal is meaningless for pixels
            outside this region.

    Returns:
    -------
    - A dictionary of boundary maps. For each dimension (e.g. if the image is 
      2D, for the first and second dimension) the dictionary contains another 
      dictionary consisting of voxel pairs. The format of the dictionary for
      each dimension is identical to the output of `get_bdry_map` (see above
      documentation). For example, output[d]['bottom']['inner'] returns the
      pixels which lie inside the excursion set at the `bottom' in dimension
      d (c.f. Fig A. in the documentation of `get_bdry_maps')
    """

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


# ============================================================================
# 
# Given a field and a threshold, the below function returns a binary image of 
# all pixels identified as neighbouring the boundary of the excursion set 
# {pixels : field at pixel = c} (this function is useful for debugging and
# visualising).
#
# ----------------------------------------------------------------------------
#
# This function takes in the following inputs:
#
# ----------------------------------------------------------------------------
#
#  - field: An image of a field  from which we wish to obtain an image of the
#           boundary of the excursion set (for example, if we want the boundary 
#           set {pixels : mean at pixel = c}, then `field' should be an image
#           of the mean.
#  - c: The value c with which to threshold the field (e.g. we are interested
#       in the set {pixels : field at pixel = c}.
#  - mask: A binary image showing pixels which should be ignored during
#          computation, or ``masked out". For example, if we are interested in
#          fMRI BOLD signal over the brain, the mask might be only the pixels
#          which lie inside the brain, as our signal is meaningless for pixels
#          outside this region.
#
# ----------------------------------------------------------------------------
#
# This function gives as outputs:
#
# ----------------------------------------------------------------------------
#
#  - A binary 1-0 image of all pixels identified as neighbouring the boundary
#    of the excursion set {pixels : field at pixel = c}
#
# ============================================================================
def get_bdry_map_combined(field, c=None, mask=None):

    # Check if the field is a string
    if isinstance(field, str):

        # Read in field
        field = read_image(field)

    # Check if the mask is a string
    if isinstance(mask, str):

        # Read in mask
        mask = read_image(mask)

    # If c is not provided, assume we have been given a binary image already
    if c is None:
        
        # Set c to 0.5
        c = 0.5

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

# ============================================================================
# 
# Given a field and a threshold, the below function returns the locations of
# all pairs of neighbouring pixels which lie inside and outside of the
# excursion set given by {pixels in field: field at pixel > c}, respectively.
#
# ----------------------------------------------------------------------------
#
# This function takes in the following inputs:
#
# ----------------------------------------------------------------------------
#
#  - bdry_maps: The dictionary of boundary maps returned by `get_bdry_maps`.
#               For each dimension (e.g. if the image is 2D, for the first and
#               second dimension) bdry_maps contains another dictionary
#               consisting of voxel pairs. The format of the dictionary for
#               each dimension is identical to the output of `get_bdry_map`
#               (see above documentation). For example:
#                   
#               bdry_maps[d]['bottom']['inner'] gives a binary image of the
#               pixels which lie inside the excursion set at the `bottom' in
#               dimension d (c.f. Fig A. in the documentation of
#               `get_bdry_maps')
#
# ----------------------------------------------------------------------------
#
# This function gives as outputs:
#
# ----------------------------------------------------------------------------
#
#  - A dictionary of voxel locations. For each dimension (e.g. if the image is 
#    2D, for the first and second dimension) the dictionary contains another 
#    dictionary consisting of locations (or indices) of voxels pairs. The
#    format of the dictionary for each dimension is designed to mirror that of 
#    the input. For example
#       
#        output[d]['bottom']['inner'] returns the locations of the pixels which
#        lie inside the excursion set at the `bottom' in dimension d (c.f.
#        Fig A. in the documentation of `get_bdry_maps')
#
# ============================================================================
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


# ============================================================================
# 
# Given a field and a dictionary of boundary locations, the below function 
# returns values of that field at the given boundary locations.
#
# ----------------------------------------------------------------------------
#
# This function takes in the following inputs:
#
# ----------------------------------------------------------------------------
#
#  - field: An image of a field  from which we wish to obtain an image of the
#           boundary of the excursion set (for example, if we want the boundary 
#           set {pixels : mean at pixel = c}, then `field' should be an image
#           of the mean.
#  - bdry_locs: A dictionary of voxel locations output by `get_bdry_locs`.
#               For each dimension (e.g. if the image is 2D, for the first and
#               second dimension) the dictionary contains another dictionary
#               consisting of locations (or indices) of voxels pairs. For
#               example:
#                    output[d]['bottom']['inner'] returns the locations of 
#                    the pixels which lie inside the excursion set at the 
#                    `bottom' in dimension d (c.f. Fig A. in the documentation 
#                    of `get_bdry_maps')
#
# ----------------------------------------------------------------------------
#
# This function gives as outputs:
#
# ----------------------------------------------------------------------------
#
#  - bdry_vals: A dictionary of voxel values. For each dimension (e.g. if the
#               image is 2D, for the first and second dimension) the dictionary
#               contains another dictionary consisting of the values of the 
#               field at the voxels pairs along the boundary. The format of 
#               the dictionary for each dimension is designed to mirror that of 
#               the bdry_locs input. For example
#                    output[d]['bottom']['inner'] returns the values of the 
#                    field at the pixels which lie inside the excursion set 
#                    at the `bottom' in dimension d (c.f. Fig A. in the
#                    documentation of `get_bdry_maps')
#
# ============================================================================
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
            inner_vals = field[(...,*bdry_locs[d][direction]['inner'])]
            outer_vals = field[(...,*bdry_locs[d][direction]['outer'])]

            # Record boundary values
            bdry_vals[d][direction]['inner'] = inner_vals
            bdry_vals[d][direction]['outer'] = outer_vals

    # Add the non-zero dimensions as an array for good measure
    bdry_vals['dims'] = bdry_locs['dims']

    # Return boundary values
    return(bdry_vals)



# ============================================================================
# 
# Given a dictionary of boundary values at pairs of voxels and a threshold,
# the below function returns interpolation weights describing where the
# threshold value lies between each pair at a sub-voxel level.
#
# ----------------------------------------------------------------------------
#
# This function takes in the following inputs:
#
# ----------------------------------------------------------------------------
#
#  - bdry_vals: A dictionary of voxel values returned by `get_bdry_vals`. For
#               each dimension (e.g. if theima ge is 2D, for the first and
#               second dimension) the dictionary contains another dictionary
#               consisting of the values of the field at the voxels pairs
#               along the boundary. The format of the dictionary for each
#               dimension is designed to mirror that of the bdry_locs input.
#               For example:
#                    output[d]['bottom']['inner'] gives the values of the 
#                    field at the pixels which lie inside the excursion set 
#                    at the `bottom' in dimension d (c.f. Fig A. in the
#                    documentation of `get_bdry_maps')
#  - c: The value c which was used to threshold the field (e.g. we are
#       interested in the set {pixels : field at pixel = c}.
#
# ----------------------------------------------------------------------------
#
# This function gives as outputs:
#
# ----------------------------------------------------------------------------
#
#   - bdry_weights: For each pair of pixels values along the boundary, p1 and 
#                   p2, this function returns a pair of weights, w1 and w2, 
#                   such that w1*p1 + w2*p2 = c. These weights represent the
#                   location of the boundary ``between the pixels'' and can be
#                   used in interpolation to get values of other fields along
#                   the boundary (c.f. Bowring (2018); Spatial confidence sets
#                   for raw effect size images).
#
#     mu
#      |                                         mu(s)=c       . straight line 
#      |                                             /     .    between pixel 1 
#      |                                            /  .          and pixel 2
#  p2 _| _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ ________ ./
#   c _| _ _ _ _ _ _ _ _ _ _ _ _ _ _ ___/ _  .  |
#      |                _______     /   .  |    |
#  p1 _| _ _ _ _ _ ____/ _ _ _ \___/.      |    |               w1 = (c-p1)/(p2-p1)
#      |         _/             . |        |    |               w2 = (p2-c)/(p2-p1)
#      |     ___/           .     |        |    |
#      |____/           .         |        |    |           
#      |____________._____________|________|____|________________
#                                 |        |    |                  S
#                              pixel 1     | pixel 2
#                                          |
#                                     Interpolated 
#                                       boundary
#                                       location
#
#        Fig B. Interpolation along the Boundary to get weights w1 and w2
#               such that w1*p1 + w2*p2 = c.
#        
# ============================================================================
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



# ============================================================================
# 
# The below function takes in a bdry_vals, bdry_weights
#
# ----------------------------------------------------------------------------
#
# This function takes in the following inputs:
#
# ----------------------------------------------------------------------------
#
#  -  field image, boundary locations
#  - boundary locations dictionary
#  - dictform
#
# ----------------------------------------------------------------------------
#
# This function gives as outputs:
#
# ----------------------------------------------------------------------------
#
#  - interpolated
#
# ============================================================================
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
            inner_vals = field[(...,*bdry_locs[d][direction]['inner'])]
            outer_vals = field[(...,*bdry_locs[d][direction]['outer'])]

            # Reshape for concatenation
            inner_vals = inner_vals.reshape(tuple(inner_vals.shape) + (1,))
            outer_vals = outer_vals.reshape(tuple(outer_vals.shape) + (1,))
            
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
        bdry_weights_inner_concat = bdry_weights_inner_concat.reshape(*bdry_weights_inner_concat.shape,1)
        bdry_weights_outer_concat = bdry_weights_outer_concat.reshape(*bdry_weights_outer_concat.shape,1)

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

