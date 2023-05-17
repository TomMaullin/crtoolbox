import os
import time
import pandas as pd
import numpy as np
from PIL import Image
import nibabel as nib

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# This file contains all miscellaneous file i/o functions used by the conf 
# sets code. These functions exist so that basic file handling does not take
# too much space in the bulk of the main code.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Tom Maullin (Last edited 10/11/2020)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_image(fname):
    """ 
    Read images from a filename.

    Parameters:
    -----------
    fname : str
        filename

    Returns:
    --------
    img : array
        image
    """

    # If the image can be loaded in with the Image package then load it in
    try:

        # Load in image
        img = Image.open(fname)

        # If the image is RGB convert it to greyscale
        if img.mode == 'RGB':

            # Convert to greyscale
            img = img.convert('L')

        # Convert to array
        img = np.array(img)

    except:

        # If the image cannot be loaded in with the Image package then try with
        # neuroimaging packages
        try:

            # Load in image
            img = nib.load(fname)

            # Convert to array
            img = np.array(img.dataobj)

        except:

            # Print error
            print("Error loading image: " + fname)

    # Return image
    return img



def read_images(fnames):
    """ 
    Read images from a list of filenames.

    Parameters:
    -----------
    fnames : list
        list of filenames

    Returns:
    --------
    imgs : array
        array of images
    """

    # Check if we have a single filename or a list of filenames
    if isinstance(fnames, str):
            
        # Convert to list
        fnames = [fnames]

    # Loop through files
    for i in range(0, len(fnames)):

        # Read in image
        img = read_image(fnames[i])

        # Check if we are on the first image
        if i == 0:

            # Record image size
            img_size = img.shape

            # Remove dimensions of length 1
            img_size = tuple([x for x in img_size if x != 1])

            # Record size of collection of images
            imgs_size = img_size + (len(fnames),)

            # Initialize array
            imgs = np.zeros(imgs_size, dtype=img.dtype)

        else:

            # Record new image size
            img_size_new = img.shape

            # Remove dimensions of length 1
            img_size_new = tuple([x for x in img_size_new if x != 1])

            # Check if image is the same size as the first image
            if img_size_new != img_size:

                # Print error
                print("Error: Images are not all the same size")

                # Exit
                exit()

        # Save image
        imgs[...,i] = img.reshape(imgs[...,i].shape)

    return imgs


def read_images_elements(fnames, indices):
    """ 
    Read images from a list of filenames one by one and retrieve the elements at
    the given indices.

    Parameters:
    -----------
    fnames : list
        List of filenames.
    indices : array
        Array of indices to retrieve.

    Returns:
    --------
    elements : array
        Array of elements at the given indices.
    """

    # Check if indices are flattened or not
    if len(indices.shape) > 1:
            
        # Flatten indices
        indices = indices.flatten()

    # Loop through files
    for i in range(0, len(fnames)):

        # Read in image
        img = read_image(fnames[i])

        # Check if we are on the first image
        if i == 0:

            # Initialize array
            elements = np.zeros((len(fnames), len(indices)))

            # Record image size
            img_size = img.shape

        else:

            # Check if image is the same size as the first image
            if img.shape != img_size:

                # Print error
                print("Error: Images are not all the same size")

                # Exit
                exit()

        # Flatten image
        img = img.flatten()

        # Save elements
        elements[i,...] = img[indices]

    return elements



def append_to_file(fname, data, remove=False):
    """
    This function takes in the data `data` and appends it to the csv file
    `fname`. If the file does not exist already it creates it. A file lock 
    system is also implemented to ensure the file is not edited by multiple
    jobs at the same time.

    Parameters:
    -----------
    - `fname`: The filename of the file we wish to append data to.
    - `data`: The data we wish to append to the file.

    Returns:
    --------
    - `None`
    """

    # Check if file is in use
    fileLocked = True
    while fileLocked:

        try:

            # Create lock file, so other jobs know we are writing to this file
            os.open(fname + ".lock", os.O_CREAT|os.O_EXCL|os.O_RDWR)
            fileLocked = False

        except FileExistsError:

        	# File is still locked
            fileLocked = True

    # Check if we are removing the file
    if remove:

        # Check if it exists
        if os.path.isfile(fname):

            # Remove it
            os.remove(fname)

    # Check whether the file exists already
    if not os.path.isfile(fname):

        # Output data to file
        pd.DataFrame(data).to_csv(fname, header=False, index=False)
        
    else:

    	# Read in data from file
    	pd.DataFrame(data).to_csv(fname, mode='a', header=False, index=False)

    # Delete lock file, so other jobs know they can now write to the
    # file
    os.remove(fname + ".lock")


def str2vec(c):
    """
    The below function takes in a string representing a vector and returns the
    vector as an array.

    Parameters:
    -----------
    - `c`: A string representing a vector.

    Returns:
    --------
    - `c`: The vector as an array.
    """

    c = str(c)
    c = c.replace("'", "")
    c = c.replace('][', '], [').replace('],[', '], [').replace('] [', '], [')
    c = c.replace('[ [', '[[').replace('] ]', ']]')
    cs = c.split(' ')
    cf = ''
    for i in range(0,len(cs)):
        cs[i]=cs[i].replace(',', '')
        cf=cf + cs[i]
        if i < (len(cs)-1):
            cf = cf + ', '
        
    return(eval(cf))


def addBlockToNifti(fname, block, blockInds,dim=None,volInd=None,aff=None,hdr=None):
    """
    The below function adds a block of voxels to a pre-existing NIFTI or creates
    a NIFTI of specified dimensions if not.

    Parameters:
    -----------
    - `fname`: An absolute path to the Nifti file.
    - `block`: The block of values to write to the NIFTI.
    - `blockInds`: The indices representing the 3D coordinates `block` should be 
                written to in the NIFTI. (Note: It is assumed if the NIFTI is
                4D we assume that the indices we want to write to in each 3D
                volume/slice are the same across all 3D volumes/slices).
    - `dim` (optional): If creating the NIFTI image for the first time, the 
                        dimensions of the NIFTI image must be specified.
    - `volInd` (optional): If we only want to write to one 3D volume/slice,
                        within a 4D file, this specifies the index of the
                        volume of interest.
    - `aff` (optional): If creating the NIFTI image for the first time, the 
                        affine of the NIFTI image must be specified.
    - `hdr` (optional): If creating the NIFTI image for the first time, the 
                        header of the NIFTI image must be specified.

    Returns:
    --------
    - `None`
    """

    # Check if file is in use
    fileLocked = True
    while fileLocked:
        try:
            # Create lock file, so other jobs know we are writing to this file
            f = os.open(fname + ".lock", os.O_CREAT|os.O_EXCL|os.O_RDWR)
            fileLocked = False
        except FileExistsError:
            fileLocked = True

    # Check volInd is correct datatype
    if volInd is not None:

        volInd = int(volInd)

    # Check whether the NIFTI exists already
    if os.path.isfile(fname):

        # Work out dim if we don't already have it
        dim = nib.Nifti1Image.from_filename(fname, mmap=False).shape

        # Work out data
        data = nib.Nifti1Image.from_filename(fname, mmap=False).get_fdata().copy()

        # Work out affine
        affine = nib.Nifti1Image.from_filename(fname, mmap=False).affine.copy()
        
    else:

        # If we know how, make the NIFTI
        if dim is not None:
            
            # Make data
            data = np.zeros(dim)

            # Make affine
            if aff is None:
                affine = np.eye(4)
            else:
                affine = aff

        else:

            # Throw an error because we don't know what to do
            raise Exception('NIFTI does not exist and dimensions not given')

    # Work out the number of output volumes inside the nifti 
    if len(dim)==3:

        # We only have one volume in this case
        n_vol = 1
        dim = np.array([dim[0],dim[1],dim[2],1])

    else:

        # The number of volumes is the last dimension
        n_vol = dim[3]

    # Seperate copy of data for outputting
    data_out = np.array(data).reshape(dim)

    # Work out the number of voxels
    n_vox = np.prod(dim[:3])

    # Reshape     
    data = data.reshape([n_vox, n_vol])

    # Add all the volumes
    if volInd is None:

        # Add block
        data[blockInds,:] = block.reshape(data[blockInds,:].shape)
        
        # Cycle through volumes, reshaping.
        for k in range(0,data.shape[1]):

            data_out[:,:,:,k] = data[:,k].reshape(int(dim[0]),
                                                  int(dim[1]),
                                                  int(dim[2]))

    # Add the one volume in the correct place
    else:

        # We're only looking at this volume
        data = data[:,volInd].reshape((n_vox,1))

        # Add block
        data[blockInds,:] = block.reshape(data[blockInds,:].shape)
        
        # Put in the volume
        data_out[:,:,:,volInd] = data[:,0].reshape(int(dim[0]),
                                                 int(dim[1]),
                                                 int(dim[2]))
    
    # Save NIFTI
    nib.save(nib.Nifti1Image(data_out, affine, header=hdr), fname)

    # Delete lock file, so other jobs know they can now write to the
    # file
    os.remove(fname + ".lock")
    os.close(f)

    del fname, data_out, affine, data, dim
