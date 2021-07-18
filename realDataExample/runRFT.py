import numpy as np
import pandas as pd
import nibabel as nib
import os
import glob
from lib.generateData import *
from lib.boundary import *
from lib.fileio import *
np.seterr(all='raise')

# ============================================================================
#
# The below function adds a block of voxels to a pre-existing NIFTI or creates
# a NIFTI of specified dimensions if not.
#
# ----------------------------------------------------------------------------
#
# This function takes the following inputs:
#
# ----------------------------------------------------------------------------
#
# - `fname`: An absolute path to the Nifti file.
# - `block`: The block of values to write to the NIFTI.
# - `blockInds`: The indices representing the 3D coordinates `block` should be 
#                written to in the NIFTI. (Note: It is assumed if the NIFTI is
#                4D we assume that the indices we want to write to in each 3D
#                volume/slice are the same across all 3D volumes/slices).
# - `dim` (optional): If creating the NIFTI image for the first time, the 
#                     dimensions of the NIFTI image must be specified.
# - `volInd` (optional): If we only want to write to one 3D volume/slice,
#                        within a 4D file, this specifies the index of the
#                        volume of interest.
# - `aff` (optional): If creating the NIFTI image for the first time, the 
#                     affine of the NIFTI image must be specified.
# - `hdr` (optional): If creating the NIFTI image for the first time, the 
#                     header of the NIFTI image must be specified.
#
# ============================================================================
def addBlockToNifti(fname, block, blockInds,dim=None,volInd=None,aff=None,hdr=None):

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

        # Load in NIFTI
        img = nib.load(fname)

        # Work out dim if we don't already have it
        dim = img.shape

        # Work out data
        data = img.get_fdata()

        # Work out affine
        affine = img.affine
        
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


    # Make NIFTI
    nifti = nib.Nifti1Image(data_out, affine, header=hdr)
    
    # Save NIFTI
    nib.save(nifti, fname)

    # Delete lock file, so other jobs know they can now write to the
    # file
    os.remove(fname + ".lock")
    os.close(f)

    del nifti, fname, data_out, affine


def runRealDat():

    # Output directory
    OutDir = '/well/nichols/users/inf852/RFT_Ttest/'

    # Get folders
    folders = glob.glob('/well/nichols/shared/HCP/Unrelated_80/??????/WM/Level2/')

    # Get IDs
    IDs = [folder.replace('/well/nichols/shared/HCP/Unrelated_80/', '').replace('/WM/Level2/', '') for folder in folders]

    # Task of interest
    taskList = ('BODY','FACE','PLACE','TOOL')

    # Slice number
    slice = 50

    # number of subjects
    nSub = len(IDs)

    # NIFTI dimensions (taken from a random nifti)
    nifdim = nib.load(os.path.join(OutDir, 'MASK_diff_BODY_' + str(IDs[0]) + '.nii')).shape

    # Loop through images getting slice
    for j, task in enumerate(taskList): 

        # Loop through number of subjects
        for i in np.arange(nSub):

            # Output filename 
            file = os.path.join(OutDir, 'COPE_diff_' + str(task) + '_' + str(IDs[i]) + '.nii')

            # Get data
            data = nib.load(file).get_data()[:,slice,:]
            data = data.reshape((1,*data.shape))

            # Concatenate files
            if i == 0:
                data_concat = data
            else:
                data_concat = np.concatenate((data_concat, data),axis=0)

            # Mask
            mask = os.path.join(OutDir, 'MASK_diff_' + str(task) + '_' + str(IDs[i]) + '.nii')
            mask = nib.load(mask).get_data()[:,slice,:]
            mask = mask.reshape((1,*mask.shape))

            # Combine masks
            if i == 0:
                mask_concat = mask
            else:
                print('before')
                print(np.any(np.isnan(mask)))
                print(np.any(np.isnan(mask_concat)))
                mask_concat = mask_concat*mask
                print('after')

            print(task)
            print(i)
        # ----------------------------------------------------------------
        # Save data
        # ----------------------------------------------------------------
        # Combine data
        if j == 0:
            datas = np.array(data_concat.reshape(1,*(data_concat.shape)))
        else:
            datas = np.concatenate((datas,data_concat.reshape(1,*(data_concat.shape))),axis=0)


    # Apply mask to data
    datas = datas*mask_concat

    # Get data shape
    data_dim = datas.shape

    # Get number of bootstraps
    nBoot = 5000

    # Get p-values
    p = [0.8,0.9,0.95]

    # Get the number of p-values we're looking at
    nPvals = len(p)

    # Get tau
    tau = 1/np.sqrt(nSub)

    # Get threshold
    c = 3.5

    # Get number of fields
    m = len(taskList)

    # Loop through number of fields
    for i in np.arange(m):
        # Get data for field i
        data = datas[i,...]
        # -------------------------------------------------------------------
        # Mean and standard deviation estimates
        # -------------------------------------------------------------------
        # Obtain mu estimate
        muHat = np.mean(data, axis=0).reshape((1,data.shape[-2],data.shape[-1]))
        # Save muHats
        if i == 0:
            muHats = np.array(muHat)
        else:
            muHats = np.concatenate((muHats,muHat),axis=0)
        # Obtain sigma
        sigma = np.std(data, axis=0).reshape(muHat.shape)
        # Save sigmas
        if i == 0:
            sigmas = np.array(sigma)
        else:
            sigmas = np.concatenate((sigmas,sigma),axis=0)

        # -------------------------------------------------------------------
        # Boundary locations for AcHati
        # -------------------------------------------------------------------
        # Get boolean maps for the boundary of AcHat
        AcHat_bdry_map = get_bdry_maps(muHat, c)

        # Get coordinates for the boundary of AcHat
        AcHat_bdry_locs = get_bdry_locs(AcHat_bdry_map)

        # Delete map as we no longer need it
        del AcHat_bdry_map

        # -------------------------------------------------------------------
        # Interpolation weights for AcHati boundary (Array version)
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for AcHati
        AcHat_bdry_vals_concat = get_bdry_values_concat(muHat, AcHat_bdry_locs)

        # Obtain the weights along the boundary for AcHati
        AcHat_bdry_weights_concat = get_bdry_weights_concat(AcHat_bdry_vals_concat, c)

        # Save boundary weights
        est_bdry_weights_concat['AcHat'+str(i+1)] = AcHat_bdry_weights_concat

        # Delete values as we no longer need them
        del AcHat_bdry_vals_concat

    # -------------------------------------------------------------------
    # Get minimum field
    # -------------------------------------------------------------------
    # This is named cap as the excursion set of the minimum field is
    # the intersection of all fields (\cap in latex)
    cap_muHat = np.amin(muHats,axis=0)

    # -------------------------------------------------------------------
    # Boundary locations for FcHat
    # -------------------------------------------------------------------
    # Get boolean map for the boundary of FcHat
    FcHat_bdry_map = get_bdry_maps(cap_muHat, c)

    # Get coordinates for the boundary of FcHat
    FcHat_bdry_locs = get_bdry_locs(FcHat_bdry_map)

    # Delete maps as we no longer need them
    del FcHat_bdry_map

    # -------------------------------------------------------------------
    # Interpolation weights for FcHat boundary (Array version)
    # -------------------------------------------------------------------
    # Obtain the values along the boundary for FcHat
    FcHat_bdry_vals_concat = get_bdry_values_concat(cap_mu, FcHat_bdry_locs)

    # Obtain the weights along the boundary for FcHat
    FcHat_bdry_weights_concat = get_bdry_weights_concat(FcHat_bdry_vals_concat, c)

    # Save boundary weights
    est_bdry_weights_concat['Fc'] = FcHat_bdry_weights_concat

    # Delete values as we no longer need them
    del FcHat_bdry_vals_concat

    # Empty dicts to store muHat
    muHat_dFcHat = {}

    # Loop through to get residuals
    for i in np.arange(m):

        # -------------------------------------------------------------------
        # Residuals along dFcHat
        # -------------------------------------------------------------------

        # Obtain residuals
        resid = (datas[i,...]-muHats[i,...])/sigmas[i,...] # NEED DATA TO BE SAVED FIRST

        # Residuals along FcHat boundary
        resid_dFcHat_concat = get_bdry_values_concat(resid, FcHat_bdry_locs)

        # Save residuals
        resids_dFcHat['field'+str(i+1)] = resid_dFcHat_concat

        # -------------------------------------------------------------------
        # MuHat along dFcHat
        # -------------------------------------------------------------------

        # Obtain MuHat along FcHat
        muHat_FcHat_bdry_concat = get_bdry_values_concat(muHats[i,...], FcHat_bdry_locs)

        # Save mu
        muHat_dFcHat['field'+str(i+1)] = muHat_FcHat_bdry_concat

    # Delete data as it is longer needed
    del datas, resid

    # -------------------------------------------------------------------
    # Boundary partitions 
    # -------------------------------------------------------------------

    # Get list of possible alphas to be considered
    alphas=list(powerset(np.arange(m)+1))

    # Locations of the dalpha boundaries
    dalphaFcHat_locs = {}

    # Loop through alphas
    for alpha in alphas:

        # Loop through possible elements of alpha
        for i in (np.arange(m)+1):

            # Get indices representing where muhati is less
            # than or equal to c if i is in alpha.
            if i in alpha:

                # Muhat
                in_dalphaFcHat_i = (muHat_dFcHat['field'+str(i)][:,1] <= c)

            # Get indices representing where mui is greater than c
            # if i is not in alpha.
            else:

                # Muhat
                in_dalphaFcHat_i = (muHat_dFcHat['field'+str(i)][:,1] > c)

            # Get a boolean index telling us whether each location in dFc
            # belongs to d^\alpha Fc.
            if i == 1:

                # Initial running product for muHat
                in_dalphaFcHat = np.array(in_dalphaFcHat_i)

            else:

                # Update running product for muhat
                in_dalphaFcHat = in_dalphaFcHat*in_dalphaFcHat_i

        # Convert boolean 1/0s into coordinates
        dalphaFcHat_loc = np.where(in_dalphaFcHat)[0]

        # Save locations
        dalphaFcHat_locs[np.array2string(alpha)] = dalphaFcHat_loc

    # -------------------------------------------------------------------
    # Get residuals and muhat along boundary partitions
    # -------------------------------------------------------------------

    # Empty dicts for residuals and muHat
    resids_dFcHat_partitioned = {}
    muHat_dFcHat_partitioned = {}

    # Loop through boundary partitions
    for alpha in alphas:

        # New empty dicts for FcHat
        resids_dalphaFcHat = {}
        muHat_dalphaFcHat = {}

        # Get dalpha locations
        dalphaFcHat_loc = dalphaFcHat_locs[np.array2string(alpha)]

        # Loop through i in alpha getting values for interpolation
        for i in alpha:

            # ------------------------------------------------------
            # Residuals and mu on estimated boundary; dFcHat
            # ------------------------------------------------------

            # Get residuals for field i along dFcHat
            residsi_dFcHat = resids_dFcHat['field'+str(i)]

            # Save residuals for field i along dalphaFc
            resids_dalphaFcHat[str(i)] = residsi_dFcHat[:,dalphaFcHat_loc,:]

            # Get mu for field i along dFcHat
            muHati_dFcHat = muHat_dFcHat['field'+str(i)]

            # Save residuals for field i along dalphaFc
            muHat_dalphaFcHat[str(i)] = muHati_dFcHat[dalphaFcHat_loc,:]

        # Save residuals for alpha
        resids_dFcHat_partitioned[np.array2string(alpha)] = resids_dalphaFcHat

        # Save mu and muhat for alpha
        muHat_dFcHat_partitioned[np.array2string(alpha)] = muHat_dalphaFcHat

    # -------------------------------------------------------------------
    # Get weights from muhat along boundary partitions for interpolation
    # -------------------------------------------------------------------

    # Empty dicts for weights
    weights_dFcHat = {}

    # Loop through boundary partitions
    for alpha in alphas:

        # New empty dicts for weights
        weights_dalphaFcHat = {}

        # Loop through i in alpha getting values for interpolation
        for i in alpha:

            # Get mui and muhati along dalpha Fc
            muHati_dalphaFcHat = muHat_dFcHat_partitioned[np.array2string(alpha)][str(i)]

            # Weights for estimated boundary
            dalphaFcHat_muHati_weights = get_bdry_weights_concat(muHati_dalphaFcHat, c)

            # Save weights
            weights_dalphaFcHat[str(i)] = dalphaFcHat_muHati_weights

        # Save weights
        weights_dFcHat[np.array2string(alpha)] = weights_dalphaFcHat

    # -------------------------------------------------------------------
    # Estimated excursion sets
    # -------------------------------------------------------------------

    # Obtain AcHat for all i
    AcHat = muHats > c

    # Obtain FcHat
    FcHat = cap_muHat > c

    # -------------------------------------------------------------------
    # Bootstrap 
    # -------------------------------------------------------------------
    # Initialize empty bootstrap stores for supremum along each boundary
    # e.g. supg_dFcHat[str(alpha)][str(i)] will give the supremum of g^i 
    #      along the boundary d\alpha F_cHat
    supg_dFcHat = {}

    # Loop through boundary partitions
    for alpha in alphas:

        # Initialize
        supg_dalphaFcHat = {}

        # Loop through i in alpha getting values for interpolation
        for i in alpha:

            # Save empty bootstrap stores
            supg_dalphaFcHat[str(i)] = np.zeros(nBoot)

        # Save empty bootstrap stores
        supg_dFcHat[np.array2string(alpha)] = supg_dalphaFcHat

    # Initalize empty bootstrap stores for the min
    min_supg_dFcHat = {}

    # Loop through boundary partitions
    for alpha in alphas:

        # Save empty bootstrap stores
        min_supg_dFcHat[np.array2string(alpha)] = np.zeros(nBoot)

    # Save empty bootstrap stores
    min_supg_dFcHat['max'] = np.zeros(nBoot)

    t1 = time.time()
    # For each bootstrap record the max of the residuals along the
    # boundary
    for b in np.arange(nBoot):

        print('Bootstrap: ', b)

        # Obtain bootstrap variables
        boot_vars = 2*np.random.randint(0,2,boot_dim)-1

        # Reshape for broadcasting purposes (extra axis refers to the fact we have
        # inner and outer boundary values in the last axes of resid_FsdGc_bdry_concat
        # and resid_FsdGcHat_bdry_concat)
        boot_vars = boot_vars.reshape((*boot_vars.shape),1)
        
        # -------------------------------------------------------------------------
        # Bootstrap residuals along dalphaFcHat
        # -------------------------------------------------------------------------
        boot_resids_dFcHat = {}

        # Loop through boundary partitions
        for alpha in alphas:

            # New empty dicts for FcHat
            boot_resids_dalphaFcHat = {}

            # New empty dict for gi along dalphaFcHat
            boot_g_dalphaFcHat = {}

            # Loop through i in alpha getting gi interpolated
            for i in alpha:

                # ------------------------------------------------------
                # Get residuals for estimated boundary
                # ------------------------------------------------------

                # Get original residuals of mui along boundary dalphaFcHat
                residsi_dalphaFcHat = resids_dFcHat_partitioned[np.array2string(alpha)][str(i)]

                # ------------------------------------------------------
                # Bootstrap residuals
                # ------------------------------------------------------

                # Multiply by rademacher variables
                boot_residsi_dalphaFcHat = boot_vars*residsi_dalphaFcHat

                # ------------------------------------------------------
                # Get gi along dalpha FcHat
                # ------------------------------------------------------

                # Sum across subjects to get the bootstrapped a values along
                # the boundary of dalphaFcHat. (Note: For some reason this is 
                # much faster if performed seperately for each of the last rows. 
                # I am still looking into why this is)
                boot_gi_dalphaFcHat = np.zeros(boot_residsi_dalphaFcHat.shape[-2:])
                boot_gi_dalphaFcHat[...,0] = np.sum(boot_residsi_dalphaFcHat[...,0], axis=0)/np.sqrt(nSub)
                boot_gi_dalphaFcHat[...,1] = np.sum(boot_residsi_dalphaFcHat[...,1], axis=0)/np.sqrt(nSub)

                # Obtain bootstrap standard deviations for muHat i along dalphaFcHat. 
                # (Note: For some reason this is much faster if performed seperately
                # for each of the last rows. I am still looking into why this is)
                boot_sigmai_dalphaFcHat = np.zeros(boot_residsi_dalphaFcHat.shape[-2:])
                boot_sigmai_dalphaFcHat[...,0] = np.std(boot_residsi_dalphaFcHat[...,0], axis=0, ddof=1)
                boot_sigmai_dalphaFcHat[...,1] = np.std(boot_residsi_dalphaFcHat[...,1], axis=0, ddof=1)

                # Divide by the boostrap standard deviation of mui
                boot_gi_dalphaFcHat = boot_gi_dalphaFcHat/boot_sigmai_dalphaFcHat

                # ------------------------------------------------------
                # Interpolate along dalpha FcHat
                # ------------------------------------------------------

                # Get weights
                dalphaFcHat_muHati_bdry_weights = weights_dFcHat[np.array2string(alpha)][str(i)]

                # Interpolation for gi along dalphaFc
                boot_gi_dalphaFcHat = get_bdry_vals_interpolated_concat(boot_gi_dalphaFcHat,dalphaFcHat_muHati_bdry_weights)

                # ------------------------------------------------------
                # Get supremum along dalphaFcHat
                # ------------------------------------------------------

                # Get the maximum of gi along dalphaFcHat if it exists
                if np.prod(boot_gi_dalphaFcHat.shape) > 0:

                    boot_g_dalphaFcHat[str(i)] = boot_gi_dalphaFcHat

            # -----------------------------------------------------------------
            # Get sup_{dalpha FcHat} |min_{alpha} g^i|
            # -----------------------------------------------------------------

            # Reset first counter
            first = True

            # Loop through i in alpha getting min(gi)
            for i in alpha:

                # Check if we have boundary values for gi
                if str(i) in boot_g_dalphaFcHat:

                    # If this is the first time initalize the min(g) array
                    if first:

                        # Initialize elementwise minimum across gi of boot_g
                        boot_ming_dalphaFcHat = boot_g_dalphaFcHat[str(i)] 

                        # We are no longer looking at the first
                        first = False

                    else:

                        # Get elementwise minimum across gi of boot_g
                        boot_ming_dalphaFcHat = np.minimum(boot_ming_dalphaFcHat,boot_g_dalphaFcHat[str(i)])

            # If we saw some boundary values (i.e. we saw at least the first i),
            # then the array was initialized
            if not first:

                # Save sup(|min(gi)|) along estimated dalpha FcHat
                boot_sup_ming_dalphaFcHat = np.max(np.abs(boot_ming_dalphaFcHat))

                # Save bootstrap result
                min_supg_dFcHat[np.array2string(alpha)][b] = boot_sup_ming_dalphaFcHat

        # -----------------------------------------------------------------
        # Get max_{P(M)} sup_{dalpha Fc} |min_{alpha} g^i| and
        # max_{P(M)} sup_{dalpha FcHat} |min_{alpha} g^i| 
        # -----------------------------------------------------------------

        # Loop through alphas taking maximum
        for alpha in alphas:

            # Check if we have boundary values for dalpha FcHat
            if np.array2string(alpha) in min_supg_dFcHat:

                # Update the minimum value we've seen
                min_supg_dFcHat['max'][b] = np.maximum(min_supg_dFcHat['max'][b],min_supg_dFcHat[np.array2string(alpha)][b])

    # -------------------------------------------------------------------
    # Obtaining a from percentiles of the max distribution
    # -------------------------------------------------------------------

    # Drop the instances where the boundary length was zero
    min_supg_dFcHat['max'] = min_supg_dFcHat['max'][min_supg_dFcHat['max']!=0]

    # If we have recorded values get their quantiles
    if (np.prod(min_supg_dFcHat['max'].shape) > 0):
        # Get the a estimates for the estimated boundary
        a_estBdry = np.percentile(min_supg_dFcHat['max'], 100*p).reshape(nPvals,1,1,1) # [pvals, 1, [1 for _ in dim>1]]
    else:
        # Set to inf by default
        a_estBdry = np.Inf*np.ones((nPvals,1,1,1))

    # Reformat them to an array form useful for boolean operation
    a_estBdry = np.concatenate((-a_estBdry,a_estBdry),axis=1)

    # -------------------------------------------------------------------
    # Get FcHat^{+/-}
    # -------------------------------------------------------------------

    # Get the statistic field which defined Achat^{+/-,i}
    g = ((muHats-c)/(sigmas*tau))

    # Take minimum over i
    stat = np.amin(g,axis=0)
    stat = stat.reshape(stat.shape[-2],stat.shape[-1])

    # Obtain FcHat^+ and FcHat^- based on a from the true boundary. This variable
    # has axes corresponding to [pvalue, plus/minus, field dimensions]
    FcHat_pm_trueBdry = stat >= a_trueBdry

    # Obtain FcHat^+ and FcHat^- based on a from the estimated boundary. This variable
    # has axes corresponding to [pvalue, plus/minus, field dimensions]
    FcHat_pm_estBdry = stat >= a_estBdry

    # -------------------------------------------------------------------
    # Get FcHat
    # -------------------------------------------------------------------
    FcHat_estBdry = stat >= 0

    # -------------------------------------------------------------------
    # Output
    # -------------------------------------------------------------------

    # Get a mask representing where the slice is
    fullmask = np.zeros(nifdim)
    fullmask[:,slice,:]=np.ones(fullmask[:,slice,:].shape)

    # Block indices
    blockInds = np.where(fullmask.reshape(np.prod(nifdim)))

    # Number of voxels per slice
    vps = nifdim[-3]*nifdim[-1]

    # Add block to nifti image for FcHat
    addBlockToNifti(os.path.join(OutDir, 'FcHat.nii'), FcHat.reshape(vps), blockInds,volInd=0,dim=NIFTIsize,aff=nifti.affine,hdr=nifti.header)

    # Loop through adding to plus/minus
    for i, pVal in enumerate(p):

        # Add block to nifti image for FcHat plus
        addBlockToNifti(os.path.join(OutDir, 'FcHatPlus_' + str(pVal) + '.nii'), FcHat_pm_estBdry[i,0,...].reshape(vps), blockInds,volInd=0,dim=NIFTIsize,aff=nifti.affine,hdr=nifti.header)

        # Add block to nifti image for FcHat minus
        addBlockToNifti(os.path.join(OutDir, 'FcHatMinus_' + str(pVal) + '.nii'), FcHat_pm_estBdry[i,1,...].reshape(vps), blockInds,volInd=0,dim=NIFTIsize,aff=nifti.affine,hdr=nifti.header)
