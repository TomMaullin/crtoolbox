import os
import numpy as np 
from crtoolbox.lib.fileio import read_image
from crtoolbox.tests.generate_ni_data import addBlockToNifti

# Powerset function, taken from:
# https://stackoverflow.com/questions/1482308/how-to-get-all-subsets-of-a-set-powerset
def powerset(s):

    # Get length of x
    x = len(s)

    # Masks is a list of `aliases' for the elements of s. I.e. for the
    # i^th element of s, the corresponding element of masks is 2^(i+1).
    masks = [1 << i for i in range(x)]
    
    # a << x is the same as a*(2^x) so we are looping over length of
    # powerset here
    for i in range(1,1 << x):

        # Yeild makes a generator obect
        # zip(masks,s) is pairing each element of s with a power of 2
        yield np.array([ss for mask, ss in zip(masks, s) if i & mask])
        

# Function to get the intersection of two masks
def intersect_masks(mask1, mask2, post_fix=None):

    # Get directory from mask1
    out_dir = os.path.dirname(mask1)
    if out_dir == '':
        out_dir = '.'

    # Read in both masks
    mask1 = 1*read_image(mask1)
    mask2 = 1*read_image(mask2)

    # Check the dimensions of the masks
    if mask1.shape != mask2.shape:
        raise ValueError("Masks must have the same dimensions.")
    
    # Compute dimensions of the masks
    D = len(mask1.shape)

    # Get the intersection of the two masks
    mask = np.logical_and(mask1, mask2)

    # If a post_fix is given, add an underscore to it
    if post_fix is not None:
        post_fix = '_' + post_fix
    else:
        post_fix = ''

    # Save the intersection mask
    if D == 3:
        
        # Write sigma to file
        addBlockToNifti(os.path.join(out_dir,"mask" + post_fix + ".nii"), mask,
                        np.arange(np.prod(mask.shape)), volInd=0,dim=mask.shape)
        
        # Return file name
        return os.path.join(out_dir,"mask" + post_fix + ".nii")
        
    elif D == 2:

        # Save to numpy array
        np.save(os.path.join(out_dir,"mask" + post_fix + ".npy"), mask)

        # Return file name
        return os.path.join(out_dir,"mask" + post_fix + ".npy")




