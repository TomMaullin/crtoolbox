import os
import numpy as np
from lib.generateData import *
from lib.boundary import *
from lib.fileio import *
import yaml
from scipy.ndimage.measurements import label

# ===========================================================================
#
# Simulation where `a' is controlled by
#
# 		a = max(sup_(dAc1\intersectAc2) |g1|, sup_(Ac1\intersectAc2d)|g2|)
#
# On two cirlces. 
#
# ---------------------------------------------------------------------------
#
# Inputs:
#
# ---------------------------------------------------------------------------
#
# - `OutDir`: Output directory.
# - `nSub`: Number of subjects.
# - `nReals`: Number of realizations.
# - `c`: threshold of interest for mu.
# - `p`: numpy array of p-values.
# - `simNo`: Number under which results will be saved
#
# ---------------------------------------------------------------------------
#
# Developers note: I needed to shorter notation but there wasn't a natural
# 				   shortening for intersection and union so I have used F and
# 				   G for shorthand for intersection and union, respectively,
#				   wherever they occur. FsdG is used in this to represent
#                  (F intersect dG) union (G intersect dF), the "sd" stands
#                  for "Symmetric difference"
#
# ===========================================================================
def SpatialSims_2mu(ipath):

    # -----------------------------------------------------------------------
    # Load in inputs 
    # -----------------------------------------------------------------------
    # Read in file
    with open(ipath, 'r') as stream:
        inputs = yaml.load(stream,Loader=yaml.FullLoader)

    # Get output directory
    OutDir = inputs['OutDir']

    # Get number of subjects
    nSub = int(inputs['nSub'])

    # Get number of simulation realizations
    nReals = int(inputs['nReals'])

    # Get Threshold
    c = np.float(inputs['c'])

    # Get p values
    p = eval(inputs['p'])

    # Get simulation number
    simNo = int(inputs['simNo'])

    # FWHM
    fwhm = str2vec(inputs['FWHM']) 

    # Get number of bootstraps
    if 'nBoot' in inputs:
        nBoot = int(inputs['nBoot'])
    else:
        nBoot = 5000 # Recommended 1e4 

    # Get number of bootstraps
    if 'tau' in inputs:
        tau = eval(inputs['tau'])
    else:
        tau = 1/np.sqrt(nSub)

    # Get specification for mu 1
    muSpec1 = inputs['mu1']

    # Reformat center and fwhm if needed
    if 'center' in muSpec1:
        muSpec1['center']=eval(muSpec1['center'])
    if 'fwhm' in muSpec1:
        muSpec1['fwhm']=eval(muSpec1['fwhm'])

    # Get specification for mu 2
    muSpec2 = inputs['mu2']

    # Reformat center and fwhm if needed
    if 'center' in muSpec2:
        muSpec2['center']=eval(muSpec2['center'])
    if 'fwhm' in muSpec2:
        muSpec2['fwhm']=eval(muSpec2['fwhm'])

    # ID for the configuration
    cfgId = inputs['cfgId']

    # -----------------------------------------------------------------------

    # Simulation directory
    simDir = os.path.join(OutDir, 'sim'+str(simNo), 'cfg' + str(cfgId))
    if not os.path.exists(simDir):
        os.mkdir(simDir)

    # Timing
    t1overall = time.time()

    # Get the number of p-values we're looking at
    nPvals = len(p)

    # Dimensions of simulated data
    data_dim = np.array([nSub, 100,100])

    # Dimensions of bootstrap variables
    boot_dim = np.array([nSub, 1]) # MARKER: COULD MAKE  np.array([batchBoot, nSub, 1])


    # Initialise array for recording the whether set violations occured, for a derived
    # from bootstrapping the true boundary and estimated boundary, respectively. The
    # result in these arrays are based on voxelwise assessment of set condition 
    # violations.
    trueBdry_success = np.zeros((nReals,nPvals))
    estBdry_success = np.zeros((nReals,nPvals))

    # Initialise array for recording the whether set violations occured, for a derived
    # from bootstrapping the true boundary and estimated boundary, respectively. The
    # result in these arrays are based on interpolation based assessment of set
    # condition violations.
    trueBdry_success_intrp = np.zeros((nReals,nPvals))
    estBdry_success_intrp = np.zeros((nReals,nPvals))

    # Loop through realizations
    for r in np.arange(nReals):

        print('r: ', r)
        # -------------------------------------------------------------------
        # Data generation
        # -------------------------------------------------------------------

        # Obtain data
        data1, mu1 = get_data(muSpec1, data_dim, fwhm)
        data2, mu2 = get_data(muSpec2, data_dim, fwhm)

        # -------------------------------------------------------------------
        # Mean and variance estimates
        # -------------------------------------------------------------------

        # Obtain mu estimate
        muHat1 = np.mean(data1, axis=0).reshape(mu1.shape)
        muHat2 = np.mean(data2, axis=0).reshape(mu1.shape)

        # Obtain sigma
        sigma1 = np.std(data1, axis=0).reshape(mu1.shape)
        sigma2 = np.std(data2, axis=0).reshape(mu1.shape)

		# -------------------------------------------------------------------
        # Boundary locations for Ac1 and Ac2
        # -------------------------------------------------------------------
        # Get boolean maps for the boundary of Ac1 and Ac2
        Ac1_bdry_maps = get_bdry_maps(mu1, c)
        Ac2_bdry_maps = get_bdry_maps(mu2, c)

        # Get coordinates for the boundary of Ac1 and Ac2
        Ac1_bdry_locs = get_bdry_locs(Ac1_bdry_maps)
        Ac2_bdry_locs = get_bdry_locs(Ac2_bdry_maps)

        # Delete maps as we no longer need them
        del Ac1_bdry_maps, Ac2_bdry_maps

        # -------------------------------------------------------------------
        # Boundary locations for Fc (= Ac1 intersect Ac2)
        # -------------------------------------------------------------------
        # Get boolean maps for the boundary of Fc
        Fc_bdry_maps = get_bdry_maps(np.minimum(mu1,mu2), c)

        # Get coordinates for the boundary of Fc
        Fc_bdry_locs = get_bdry_locs(Fc_bdry_maps)

        # Delete maps as we no longer need them
        del Fc_bdry_maps

        # -------------------------------------------------------------------
        # Boundary locations for AcHat1 and AcHat2
        # -------------------------------------------------------------------
        # Get boolean maps for the boundary of AcHat1 and AcHat2
        AcHat1_bdry_maps = get_bdry_maps(muHat1, c)
        AcHat2_bdry_maps = get_bdry_maps(muHat2, c)

        # Get coordinates for the boundary of AcHat1 and AcHat2
        AcHat1_bdry_locs = get_bdry_locs(AcHat1_bdry_maps)
        AcHat2_bdry_locs = get_bdry_locs(AcHat2_bdry_maps)

        # Delete maps as we no longer need them
        del AcHat1_bdry_maps, AcHat2_bdry_maps 

        # -------------------------------------------------------------------
        # Boundary locations for FcHat(= AcHat1 intersect AcHat2)
        # -------------------------------------------------------------------
        # Get boolean maps for the boundary of FcHat
        FcHat_bdry_maps = get_bdry_maps(np.minimum(muHat1,muHat2), c)

        # Get coordinates for the boundary of FcHat
        FcHat_bdry_locs = get_bdry_locs(FcHat_bdry_maps)

        # Delete maps as we no longer need them
        del FcHat_bdry_maps

        # -------------------------------------------------------------------
        # Interpolation weights for Ac1 and Ac2 boundary (Dict version)
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for Ac1 and Ac2
        Ac1_bdry_vals = get_bdry_values(mu1, Ac1_bdry_locs)
        Ac2_bdry_vals = get_bdry_values(mu2, Ac2_bdry_locs)

        # Obtain the weights along the boundary for Ac1 and Ac2
        Ac1_bdry_weights = get_bdry_weights(Ac1_bdry_vals, c)
        Ac2_bdry_weights = get_bdry_weights(Ac2_bdry_vals, c)

        # Delete values as we no longer need them
        del Ac1_bdry_vals, Ac2_bdry_vals

        # -------------------------------------------------------------------
        # Interpolation weights for Fc boundary (Dict version)
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for Fc
        Fc_bdry_vals = get_bdry_values(np.minimum(mu1,mu2), Fc_bdry_locs)

        # Obtain the weights along the boundary for Fc
        Fc_bdry_weights = get_bdry_weights(Fc_bdry_vals, c)

        # Delete values as we no longer need them
        del Fc_bdry_vals

        # -------------------------------------------------------------------
        # Interpolation weights for Ac1 and Ac2 boundary (Array version)
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for Ac1 and Ac2
        Ac1_bdry_vals_concat = get_bdry_values_concat(mu1, Ac1_bdry_locs)
        Ac2_bdry_vals_concat = get_bdry_values_concat(mu2, Ac2_bdry_locs)

        # Obtain the weights along the boundary for Ac1 and Ac2
        Ac1_bdry_weights_concat = get_bdry_weights_concat(Ac1_bdry_vals_concat, c)
        Ac2_bdry_weights_concat = get_bdry_weights_concat(Ac2_bdry_vals_concat, c)

        # Delete values as we no longer need them
        del Ac1_bdry_vals_concat, Ac2_bdry_vals_concat

        # -------------------------------------------------------------------
        # Interpolation weights for Fc boundary (Array version)
        # -------------------------------------------------------------------
		# Obtain the values along the boundary for Fc
        Fc_bdry_vals_concat = get_bdry_values_concat(np.minimum(mu1,mu2), Fc_bdry_locs)

        # Obtain the weights along the boundary for Fc
        Fc_bdry_weights_concat = get_bdry_weights_concat(Fc_bdry_vals_concat, c)
 
 		# Delete values as we no longer need them
        del Fc_bdry_vals_concat

        # -------------------------------------------------------------------
        # Interpolation weights for AcHat1 and AcHat2 boundary (Array
        # version)
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for AcHat1 and AcHat2
        AcHat1_bdry_vals_concat = get_bdry_values_concat(muHat1, AcHat1_bdry_locs)
        AcHat2_bdry_vals_concat = get_bdry_values_concat(muHat2, AcHat2_bdry_locs)

        # Obtain the weights along the boundary for AcHat1 and AcHat2
        AcHat1_bdry_weights_concat = get_bdry_weights_concat(AcHat1_bdry_vals_concat, c)
        AcHat2_bdry_weights_concat = get_bdry_weights_concat(AcHat2_bdry_vals_concat, c)

        # Delete values as we no longer need them
        del AcHat1_bdry_vals_concat, AcHat2_bdry_vals_concat

        # -------------------------------------------------------------------
        # Interpolation weights for FcHat boundary (Array version)
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for FcHat
        FcHat_bdry_vals_concat = get_bdry_values_concat(np.minimum(muHat1,muHat2), FcHat_bdry_locs)

        # Obtain the weights along the boundary for FcHat
        FcHat_bdry_weights_concat = get_bdry_weights_concat(FcHat_bdry_vals_concat, c)

        # Delete values as we no longer need them
        del FcHat_bdry_vals_concat

        # -------------------------------------------------------------------
        # Residuals along dFcHat
        # -------------------------------------------------------------------

        # Obtain residuals
        resid1 = (data1-muHat1)/sigma1
        resid2 = (data2-muHat2)/sigma2

		# Residuals along Fc boundary
        resid1_Fc_bdry_concat = get_bdry_values_concat(resid1, Fc_bdry_locs)
        resid2_Fc_bdry_concat = get_bdry_values_concat(resid2, Fc_bdry_locs)

        # Residuals along FcHat boundary
        resid1_FcHat_bdry_concat = get_bdry_values_concat(resid1, FcHat_bdry_locs)
        resid2_FcHat_bdry_concat = get_bdry_values_concat(resid2, FcHat_bdry_locs)

        # Delete residuals as they are no longer needed
        #del data1, data2

        # -------------------------------------------------------------------
        # Mu along FsdGcHat and MuHat along FsdGcHat
        # -------------------------------------------------------------------

        # Obtain Mu along Fc
        mu1_Fc_bdry_concat = get_bdry_values_concat(mu1, Fc_bdry_locs)
        mu2_Fc_bdry_concat = get_bdry_values_concat(mu2, Fc_bdry_locs)

        # -------------------------------------------------------------------
        # Images
        # -------------------------------------------------------------------
        if r==1 and inputs['figGen']:
            
            # Make image directory
            figDir = os.path.join(OutDir, 'sim'+str(simNo), 'Figures')
            if not os.path.exists(figDir):
                os.mkdir(figDir)

            # AcHat1
            AcHat1_im = muHat1>c
            plt.figure(0)
            plt.imshow(1*AcHat1_im[0,:,:])
            plt.savefig(os.path.join(figDir, 'AcHat1_cfg'+str(cfgId)+'.png'))

            # dAcHat1
            dAcHat1_im = get_bdry_map_combined(muHat1, c)
            plt.figure(1)
            plt.imshow(1*dAcHat1_im[0,:,:])
            plt.savefig(os.path.join(figDir, 'dAcHat1_cfg'+str(cfgId)+'.png'))

            # AcHat2
            AcHat2_im = muHat2>c
            plt.figure(2)
            plt.imshow(1*AcHat2_im[0,:,:])
            plt.savefig(os.path.join(figDir, 'AcHat2_cfg'+str(cfgId)+'.png'))

            # dAcHat2
            dAcHat2_im = get_bdry_map_combined(muHat2, c)
            plt.figure(3)
            plt.imshow(1*dAcHat2_im[0,:,:])
            plt.savefig(os.path.join(figDir, 'dAcHat2_cfg'+str(cfgId)+'.png'))

            # AcHat1 \cup AcHat2
            AcHat1cupAcHat2_im = np.maximum(muHat1,muHat2)>c
            plt.figure(4)
            plt.imshow(1*AcHat1cupAcHat2_im[0,:,:])
            plt.savefig(os.path.join(figDir, 'GcHat_cfg'+str(cfgId)+'.png'))

            # d(AcHat1 \cup AcHat2)
            dAcHat1cupAcHat2_im = get_bdry_map_combined(np.maximum(muHat1,muHat2),c)
            plt.figure(5)
            plt.imshow(1*dAcHat1cupAcHat2_im[0,:,:])
            plt.savefig(os.path.join(figDir, 'dGcHat_cfg'+str(cfgId)+'.png'))

            # AcHat1 \cap AcHat2
            AcHat1capAcHat2_im = np.minimum(muHat1,muHat2)>c
            plt.figure(6)
            plt.imshow(1*AcHat1capAcHat2_im[0,:,:])
            plt.savefig(os.path.join(figDir, 'FcHat_cfg'+str(cfgId)+'.png'))

            # d(AcHat1 \cap AcHat2)
            dAcHat1capAcHat2_im = get_bdry_map_combined(np.minimum(muHat1,muHat2),c)
            plt.figure(7)
            plt.imshow(1*dAcHat1capAcHat2_im[0,:,:])
            plt.savefig(os.path.join(figDir, 'dFcHat_cfg'+str(cfgId)+'.png'))

            # dF12cHat = dAcHat1 \cap dAcHat2
            dAcHat1capdAcHat2_im = dAcHat1capAcHat2_im*dAcHat1cupAcHat2_im
            plt.figure(8)
            plt.imshow(1*dAcHat1capdAcHat2_im[0,:,:])
            plt.savefig(os.path.join(figDir, 'dF12cHat_cfg'+str(cfgId)+'.png'))

            # dF1cHat = dAcHat1 \cap AcHat2
            rhs_im = dAcHat1_im*AcHat2_im
            plt.figure(10)
            plt.imshow(1*rhs_im[0,:,:])
            plt.savefig(os.path.join(figDir, 'dF1cHat_cfg'+str(cfgId)+'.png'))

            # dF2cHat = dAcHat2 \cap AcHat1
            lhs_im = dAcHat2_im*AcHat1_im
            plt.figure(9)
            plt.imshow(1*lhs_im[0,:,:])
            plt.savefig(os.path.join(figDir, 'dF2cHat_cfg'+str(cfgId)+'.png'))

            # muHat1
            plt.figure(11)
            plt.imshow(muHat1[0,:,:])
            plt.colorbar()
            plt.savefig(os.path.join(figDir, 'muHat1_cfg'+str(cfgId)+'.png'))

            # muHat2
            plt.figure(12)
            plt.imshow(muHat2[0,:,:])
            plt.colorbar()
            plt.savefig(os.path.join(figDir, 'muHat2_cfg'+str(cfgId)+'.png'))

            # data1
            plt.figure(13)
            plt.imshow(data1[10,:,:])
            plt.colorbar()
            plt.savefig(os.path.join(figDir, 'data1_cfg'+str(cfgId)+'.png'))

            # data2
            plt.figure(14)
            plt.imshow(data2[10,:,:])
            plt.colorbar()
            plt.savefig(os.path.join(figDir, 'data2_cfg'+str(cfgId)+'.png'))

            # max(muHat1,muHat2)
            plt.figure(15)
            plt.imshow(np.maximum(muHat1[0,:,:],muHat2[0,:,:]))
            plt.colorbar()
            plt.savefig(os.path.join(figDir, 'maxMuHat_cfg'+str(cfgId)+'.png'))

            # min(muHat1,muHat2)
            plt.figure(16)
            plt.imshow(np.minimum(muHat1[0,:,:],muHat2[0,:,:]))
            plt.colorbar()
            plt.savefig(os.path.join(figDir, 'minMuHat_cfg'+str(cfgId)+'.png'))


        print('shapes here : ', mu1_Fc_bdry_concat.shape,resid1_Fc_bdry_concat.shape)

        # Get locations where outer mu1 and mu2 are greater than c
        mu1_Fc_bdry_gcloc = np.where(mu1_Fc_bdry_concat[0,:,1]>c)[0]
        mu2_Fc_bdry_gcloc = np.where(mu2_Fc_bdry_concat[0,:,1]>c)[0]

        print(resid1_Fc_bdry_concat[:,mu1_Fc_bdry_gcloc,:].shape)
        print(resid2_Fc_bdry_concat[:,mu2_Fc_bdry_gcloc,:].shape)

        # Obtain MuHat along Fc
        muHat1_FcHat_bdry_concat = get_bdry_values_concat(muHat1, FcHat_bdry_locs)
        muHat2_FcHat_bdry_concat = get_bdry_values_concat(muHat2, FcHat_bdry_locs)

        print('shapes here : ', muHat1_FcHat_bdry_concat.shape,resid1_FcHat_bdry_concat.shape)

        # Get locations where outer muHat1 and muHat2 are greater than c
        muHat1_FcHat_bdry_gcloc = np.where(muHat1_FcHat_bdry_concat[0,:,1]>c)[0]
        muHat2_FcHat_bdry_gcloc = np.where(muHat2_FcHat_bdry_concat[0,:,1]>c)[0]

        print(resid1_FcHat_bdry_concat[:,muHat1_FcHat_bdry_gcloc,:].shape)
        print(resid2_FcHat_bdry_concat[:,muHat2_FcHat_bdry_gcloc,:].shape)

        # -------------------------------------------------------------------
        # Residuals along FsdG and FsdGHat boundaries (Array version)
        # -------------------------------------------------------------------
        # In this simulation we are bootstrapping the residuals for field 1
        # along dAc1 intersect Ac2 i.e. along dF where mu2 > c. And vice versa
        # for field 2.
        resid1_FsdGc1_bdry_concat = resid1_Fc_bdry_concat[:,mu2_Fc_bdry_gcloc,:]
        resid2_FsdGc2_bdry_concat = resid2_Fc_bdry_concat[:,mu1_Fc_bdry_gcloc,:]

        # In this simulation we are bootstrapping the residuals for field 1
        # along dAcHat1 intersect AcHat2 i.e. along dF where muHat2 > c. And 
        # vice versa for field 2.
        resid1_FsdGcHat1_bdry_concat = resid1_FcHat_bdry_concat[:,muHat2_FcHat_bdry_gcloc,:]
        resid2_FsdGcHat2_bdry_concat = resid2_FcHat_bdry_concat[:,muHat1_FcHat_bdry_gcloc,:]

        # -------------------------------------------------------------------
        # Mu and MuHat along FsdG and FsdGHat boundaries (Array version)
        # -------------------------------------------------------------------
        # In this simulation we are bootstrapping the residuals for field 1
        # along dAc1 intersect Ac2 i.e. along dF where mu2 > c. And vice versa
        # for field 2.
        mu1_FsdGc1_bdry_concat = mu1_Fc_bdry_concat[:,mu2_Fc_bdry_gcloc,:]
        mu2_FsdGc2_bdry_concat = mu2_Fc_bdry_concat[:,mu1_Fc_bdry_gcloc,:]

        # In this simulation we are bootstrapping the residuals for field 1
        # along dAcHat1 intersect AcHat2 i.e. along dF where muHat2 > c. And 
        # vice versa for field 2.
        muHat1_FsdGcHat1_bdry_concat = muHat1_FcHat_bdry_concat[:,muHat2_FcHat_bdry_gcloc,:]
        muHat2_FsdGcHat2_bdry_concat = muHat2_FcHat_bdry_concat[:,muHat1_FcHat_bdry_gcloc,:]

        # -------------------------------------------------------------------
        # Interpolation weights FsdG and FsdGHat boundaries (Array version)
        # -------------------------------------------------------------------
        # Obtain the weights along the boundary for FsdGc
        FsdGc1_bdry_weights_concat = get_bdry_weights_concat(mu1_FsdGc1_bdry_concat, c)
        FsdGc2_bdry_weights_concat = get_bdry_weights_concat(mu2_FsdGc2_bdry_concat, c)

        # Obtain the weights along the boundary for FsdGcHat
        FsdGcHat1_bdry_weights_concat = get_bdry_weights_concat(muHat1_FsdGcHat1_bdry_concat, c)
        FsdGcHat2_bdry_weights_concat = get_bdry_weights_concat(muHat2_FsdGcHat2_bdry_concat, c)

        # -------------------------------------------------------------------
        # True and estimated excursion sets
        # -------------------------------------------------------------------
        # Obtain Ac1 and Ac2
        Ac1 = mu1 > c
        Ac2 = mu2 > c

        # Obtain Fc
        Fc = np.minimum(mu1,mu2) > c

        # Obtain AcHat2 and AcHat2 
        AcHat1 = muHat1 > c
        AcHat2 = muHat2 > c

        # Obtain FcHat
        FcHat = np.minimum(muHat1,muHat2) > c

        # -------------------------------------------------------------------
        # Bootstrap 
        # -------------------------------------------------------------------
        # Initialize empty bootstrap stores
        max_g1_FsdGc1 = np.zeros(nBoot)
        max_g2_FsdGc2 = np.zeros(nBoot)
        max_g1_FsdGcHat1 = np.zeros(nBoot)
        max_g2_FsdGcHat2 = np.zeros(nBoot)

        t1 = time.time()
        # For each bootstrap record the max of the residuals along the
        # boundary
        for b in np.arange(nBoot):

            # Obtain bootstrap variables
            boot_vars = 2*np.random.randint(0,2,boot_dim)-1

            # Reshape for broadcasting purposes (extra axis refers to the fact we have
            # inner and outer boundary values in the last axes of resid_FsdGc_bdry_concat
            # and resid_FsdGcHat_bdry_concat)
            boot_vars = boot_vars.reshape((*boot_vars.shape),1)

            # Bootstrap residuals along FsdGc1 and FsdGc2
            boot_resid1_FsdGc1_bdry_concat = boot_vars*resid1_FsdGc1_bdry_concat
            boot_resid2_FsdGc2_bdry_concat = boot_vars*resid2_FsdGc2_bdry_concat

            # Bootstrap residuals along FsdGcHat1 and FsdGcHat2
            boot_resid1_FsdGcHat1_bdry_concat = boot_vars*resid1_FsdGcHat1_bdry_concat
            boot_resid2_FsdGcHat2_bdry_concat = boot_vars*resid2_FsdGcHat2_bdry_concat

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of FsdGc1. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g1_FsdGc1_bdry_concat = np.zeros(boot_resid1_FsdGc1_bdry_concat.shape[-2:])
            boot_g1_FsdGc1_bdry_concat[...,0] = np.sum(boot_resid1_FsdGc1_bdry_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g1_FsdGc1_bdry_concat[...,1] = np.sum(boot_resid1_FsdGc1_bdry_concat[...,1], axis=0)/np.sqrt(nSub)

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of FsdGc2. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g2_FsdGc2_bdry_concat = np.zeros(boot_resid2_FsdGc2_bdry_concat.shape[-2:])
            boot_g2_FsdGc2_bdry_concat[...,0] = np.sum(boot_resid2_FsdGc2_bdry_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g2_FsdGc2_bdry_concat[...,1] = np.sum(boot_resid2_FsdGc2_bdry_concat[...,1], axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along FsdGc1. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma1_boot_FsdGc1_concat = np.zeros(boot_resid1_FsdGc1_bdry_concat.shape[-2:])
            sigma1_boot_FsdGc1_concat[...,0] = np.std(boot_resid1_FsdGc1_bdry_concat[...,0], axis=0, ddof=1)
            sigma1_boot_FsdGc1_concat[...,1] = np.std(boot_resid1_FsdGc1_bdry_concat[...,1], axis=0, ddof=1)

            # Obtain bootstrap standard deviations along FsdGc2. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma2_boot_FsdGc2_concat = np.zeros(boot_resid2_FsdGc2_bdry_concat.shape[-2:])
            sigma2_boot_FsdGc2_concat[...,0] = np.std(boot_resid2_FsdGc2_bdry_concat[...,0], axis=0, ddof=1)
            sigma2_boot_FsdGc2_concat[...,1] = np.std(boot_resid2_FsdGc2_bdry_concat[...,1], axis=0, ddof=1)

            # Divide by the boostrap standard deviation on FsdGc1
            boot_g1_FsdGc1_bdry_concat = boot_g1_FsdGc1_bdry_concat/sigma1_boot_FsdGc1_concat

            # Divide by the boostrap standard deviation on FsdGc2
            boot_g2_FsdGc2_bdry_concat = boot_g2_FsdGc2_bdry_concat/sigma2_boot_FsdGc2_concat

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of FsdGcHat1. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g1_FsdGcHat1_bdry_concat = np.zeros(boot_resid1_FsdGcHat1_bdry_concat.shape[-2:])
            boot_g1_FsdGcHat1_bdry_concat[...,0] = np.sum(boot_resid1_FsdGcHat1_bdry_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g1_FsdGcHat1_bdry_concat[...,1] = np.sum(boot_resid1_FsdGcHat1_bdry_concat[...,1], axis=0)/np.sqrt(nSub)

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of FsdGcHat2. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g2_FsdGcHat2_bdry_concat = np.zeros(boot_resid2_FsdGcHat2_bdry_concat.shape[-2:])
            boot_g2_FsdGcHat2_bdry_concat[...,0] = np.sum(boot_resid2_FsdGcHat2_bdry_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g2_FsdGcHat2_bdry_concat[...,1] = np.sum(boot_resid2_FsdGcHat2_bdry_concat[...,1], axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along FsdGcHat1. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma1_boot_FsdGcHat1_concat = np.zeros(boot_resid1_FsdGcHat1_bdry_concat.shape[-2:])
            sigma1_boot_FsdGcHat1_concat[...,0] = np.std(boot_resid1_FsdGcHat1_bdry_concat[...,0], axis=0, ddof=1)
            sigma1_boot_FsdGcHat1_concat[...,1] = np.std(boot_resid1_FsdGcHat1_bdry_concat[...,1], axis=0, ddof=1)

            # Obtain bootstrap standard deviations along FsdGcHat2. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma2_boot_FsdGcHat2_concat = np.zeros(boot_resid2_FsdGcHat2_bdry_concat.shape[-2:])
            sigma2_boot_FsdGcHat2_concat[...,0] = np.std(boot_resid2_FsdGcHat2_bdry_concat[...,0], axis=0, ddof=1)
            sigma2_boot_FsdGcHat2_concat[...,1] = np.std(boot_resid2_FsdGcHat2_bdry_concat[...,1], axis=0, ddof=1)

            # Divide by the boostrap standard deviation on FsdGcHat1
            boot_g1_FsdGcHat1_bdry_concat = boot_g1_FsdGcHat1_bdry_concat/sigma1_boot_FsdGcHat1_concat

            # Divide by the boostrap standard deviation on FsdGcHat2
            boot_g2_FsdGcHat2_bdry_concat = boot_g2_FsdGcHat2_bdry_concat/sigma2_boot_FsdGcHat2_concat

            # Interpolation for FsdGc1 and FsdGc2 boundary
            boot_g1_FsdGc1_bdry_concat = get_bdry_vals_interpolated_concat(boot_g1_FsdGc1_bdry_concat,FsdGc1_bdry_weights_concat)
            boot_g2_FsdGc2_bdry_concat = get_bdry_vals_interpolated_concat(boot_g2_FsdGc2_bdry_concat,FsdGc2_bdry_weights_concat)

            # Interpolation for FsdGcHat1 and FsdGcHat2 boundary
            boot_g1_FsdGcHat1_bdry_concat = get_bdry_vals_interpolated_concat(boot_g1_FsdGcHat1_bdry_concat,FsdGcHat1_bdry_weights_concat)
            boot_g2_FsdGcHat2_bdry_concat = get_bdry_vals_interpolated_concat(boot_g2_FsdGcHat2_bdry_concat,FsdGcHat2_bdry_weights_concat)

            # Get maximum along FsdGc1 and FsdGc2 boudary
            try:
                max_g1_FsdGc1[b] = np.max(np.abs(boot_g1_FsdGc1_bdry_concat)) 
                max_g2_FsdGc2[b] = np.max(np.abs(boot_g2_FsdGc2_bdry_concat)) 

                # Get maximum along FsdGcHat1 and FsdGcHat2 boudary
                max_g1_FsdGcHat1[b] = np.max(np.abs(boot_g1_FsdGcHat1_bdry_concat))
                max_g2_FsdGcHat2[b] = np.max(np.abs(boot_g2_FsdGcHat2_bdry_concat)) 

            except:
            
                # Make image directory
                fEDir = os.path.join(OutDir, 'sim'+str(simNo), 'FiguresErrored')
                if not os.path.exists(fEDir):
                    os.mkdir(fEDir)

                # AcHat1
                AcHat1_im = muHat1>c
                plt.figure(0)
                plt.imshow(1*AcHat1_im[0,:,:])
                plt.savefig(os.path.join(fEDir, 'AcHat1_cfg'+str(cfgId)+'.png'))

                # dAcHat1
                dAcHat1_im = get_bdry_map_combined(muHat1, c)
                plt.figure(1)
                plt.imshow(1*dAcHat1_im[0,:,:])
                plt.savefig(os.path.join(fEDir, 'dAcHat1_cfg'+str(cfgId)+'.png'))

                # AcHat2
                AcHat2_im = muHat2>c
                plt.figure(2)
                plt.imshow(1*AcHat2_im[0,:,:])
                plt.savefig(os.path.join(fEDir, 'AcHat2_cfg'+str(cfgId)+'.png'))

                # dAcHat2
                dAcHat2_im = get_bdry_map_combined(muHat2, c)
                plt.figure(3)
                plt.imshow(1*dAcHat2_im[0,:,:])
                plt.savefig(os.path.join(fEDir, 'dAcHat2_cfg'+str(cfgId)+'.png'))

                # AcHat1 \cup AcHat2
                AcHat1cupAcHat2_im = np.maximum(muHat1,muHat2)>c
                plt.figure(4)
                plt.imshow(1*AcHat1cupAcHat2_im[0,:,:])
                plt.savefig(os.path.join(fEDir, 'GcHat_cfg'+str(cfgId)+'.png'))

                # d(AcHat1 \cup AcHat2)
                dAcHat1cupAcHat2_im = get_bdry_map_combined(np.maximum(muHat1,muHat2),c)
                plt.figure(5)
                plt.imshow(1*dAcHat1cupAcHat2_im[0,:,:])
                plt.savefig(os.path.join(fEDir, 'dGcHat_cfg'+str(cfgId)+'.png'))

                # AcHat1 \cap AcHat2
                AcHat1capAcHat2_im = np.minimum(muHat1,muHat2)>c
                plt.figure(6)
                plt.imshow(1*AcHat1capAcHat2_im[0,:,:])
                plt.savefig(os.path.join(fEDir, 'FcHat_cfg'+str(cfgId)+'.png'))

                # d(AcHat1 \cap AcHat2)
                dAcHat1capAcHat2_im = get_bdry_map_combined(np.minimum(muHat1,muHat2),c)
                plt.figure(7)
                plt.imshow(1*dAcHat1capAcHat2_im[0,:,:])
                plt.savefig(os.path.join(fEDir, 'dFcHat_cfg'+str(cfgId)+'.png'))

                # dF12cHat = dAcHat1 \cap dAcHat2
                dAcHat1capdAcHat2_im = dAcHat1capAcHat2_im*dAcHat1cupAcHat2_im
                plt.figure(8)
                plt.imshow(1*dAcHat1capdAcHat2_im[0,:,:])
                plt.savefig(os.path.join(fEDir, 'dF12cHat_cfg'+str(cfgId)+'.png'))

                # dF1cHat = dAcHat1 \cap AcHat2
                rhs_im = dAcHat1_im*AcHat2_im
                plt.figure(10)
                plt.imshow(1*rhs_im[0,:,:])
                plt.savefig(os.path.join(fEDir, 'dF1cHat_cfg'+str(cfgId)+'.png'))

                # dF2cHat = dAcHat2 \cap AcHat1
                lhs_im = dAcHat2_im*AcHat1_im
                plt.figure(9)
                plt.imshow(1*lhs_im[0,:,:])
                plt.savefig(os.path.join(fEDir, 'dF2cHat_cfg'+str(cfgId)+'.png'))

                # muHat1
                plt.figure(11)
                plt.imshow(muHat1[0,:,:])
                plt.colorbar()
                plt.savefig(os.path.join(fEDir, 'muHat1_cfg'+str(cfgId)+'.png'))

                # muHat2
                plt.figure(12)
                plt.imshow(muHat2[0,:,:])
                plt.colorbar()
                plt.savefig(os.path.join(fEDir, 'muHat2_cfg'+str(cfgId)+'.png'))

                # data1
                plt.figure(13)
                plt.imshow(data1[10,:,:])
                plt.colorbar()
                plt.savefig(os.path.join(fEDir, 'data1_cfg'+str(cfgId)+'.png'))

                # data2
                plt.figure(14)
                plt.imshow(data2[10,:,:])
                plt.colorbar()
                plt.savefig(os.path.join(fEDir, 'data2_cfg'+str(cfgId)+'.png'))

                # max(muHat1,muHat2)
                plt.figure(15)
                plt.imshow(np.maximum(muHat1[0,:,:],muHat2[0,:,:]))
                plt.colorbar()
                plt.savefig(os.path.join(fEDir, 'maxMuHat_cfg'+str(cfgId)+'.png'))

                # min(muHat1,muHat2)
                plt.figure(16)
                plt.imshow(np.minimum(muHat1[0,:,:],muHat2[0,:,:]))
                plt.colorbar()
                plt.savefig(os.path.join(fEDir, 'minMuHat_cfg'+str(cfgId)+'.png'))


            # -------------------------------------------------------------------------------------

        # Get the maximum needed for the true boundary
        max_g_FsdGc = np.maximum(max_g1_FsdGc1,max_g2_FsdGc2)

        # Get the maximum needed for the estimated boundary
        max_g_FsdGcHat = np.maximum(max_g1_FsdGcHat1,max_g2_FsdGcHat2)

        t2 = time.time()
        print('Bootstrap time: ', t2-t1)

        # -------------------------------------------------------------------
        # Obtaining a from percentiles of the max distribution
        # -------------------------------------------------------------------

        # Get the a estimates for the true boundaries and estimated boundaries
        a_trueBdry = np.percentile(max_g_FsdGc, 100*p).reshape(nPvals,1,1,1)
        a_estBdry = np.percentile(max_g_FsdGcHat, 100*p).reshape(nPvals,1,1,1) # [pvals, 1, [1 for _ in dim>1]]

        # Reformat them to an array form useful for boolean operation
        a_trueBdry = np.concatenate((-a_trueBdry,a_trueBdry),axis=1)
        a_estBdry = np.concatenate((-a_estBdry,a_estBdry),axis=1)

        # Get the statistic field which defined Achat^{+/-,1/2}
        stat1 = ((muHat1-c)/(sigma1*tau)).reshape(1,(*muHat1.shape))
        stat2 = ((muHat2-c)/(sigma2*tau)).reshape(1,(*muHat2.shape))
        
        # Minimum for intersection
        stat = np.minimum(stat1,stat2)
        stat = stat.reshape(stat.shape[-2],stat.shape[-1])

        # Obtain FcHat^+ and FcHat^- based on a from the true boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        FcHat_pm_trueBdry = stat >= a_trueBdry

        # Obtain FcHat^+ and FcHat^- based on a from the estimated boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        FcHat_pm_estBdry = stat >= a_estBdry

        # -------------------------------------------------------------------
        # Some set logic to work out violations
        # -------------------------------------------------------------------

        # Obtain FcHat^+\FcHat based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        FcHatp_sub_Fc_trueBdry = FcHat_pm_trueBdry[:,1,...] & ~Fc[...]

        # Obtain FcHat^+\FcHat based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        FcHatp_sub_Fc_estBdry = FcHat_pm_estBdry[:,1,...] & ~Fc[...]

        # Obtain Fc\FcHat^- based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        Fc_sub_FcHatm_trueBdry = Fc[...] & ~FcHat_pm_trueBdry[:,0,...]

        # Obtain Fc\FcHat^- based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        Fc_sub_FcHatm_estBdry = Fc[...] & ~FcHat_pm_estBdry[:,0,...]

        #print(Fc_sub_FcHatm_estBdry.shape, Fc.shape, FcHat_pm_estBdry.shape)

        # -------------------------------------------------------------------
        # Check whether there were any boundary violations using voxelwise
        # set logic (checking if voxels existed in one set but not another, 
        # etc)
        # -------------------------------------------------------------------
        # Record if we saw a violation in the true boundary based sets
        trueBdry_success[r,:] = 1-(np.any(FcHatp_sub_Fc_trueBdry,axis=(1,2)) | np.any(Fc_sub_FcHatm_trueBdry,axis=(1,2))) # MARKER: AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success[r,:] = 1-(np.any(FcHatp_sub_Fc_estBdry,axis=(1,2)) | np.any(Fc_sub_FcHatm_estBdry,axis=(1,2)))# MARKER: AXES WONT WORK FOR 3D ATM

        # -------------------------------------------------------------------
        # Get stat along the Fc boundary
        # -------------------------------------------------------------------
        # Get the values along the outer and inner boundaries
        stat_FcBdry = get_bdry_values(stat, Fc_bdry_locs)

        # Interpolate to get the values along the true boundary
        stat_FcBdry = get_bdry_vals_interpolated(stat_FcBdry, Fc_bdry_weights)

        # -------------------------------------------------------------------
        # Check whether there were any boundary violations using interpolated
        # boundary values (checking if voxels had values corresponding to no
        # violations, etc)
        # -------------------------------------------------------------------

        # Perform lower check on stat map using thresholds based on the
        # estimated boundary
        bdry_lowerCheck_estBdry = stat_FcBdry >= a_estBdry[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # estimated boundary
        bdry_upperCheck_estBdry = stat_FcBdry <= a_estBdry[:,1,:,0]

        # Perform lower check on stat map using thresholds based on the
        # true boundary
        bdry_lowerCheck_trueBdry = stat_FcBdry >= a_trueBdry[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # true boundary
        bdry_upperCheck_trueBdry = stat_FcBdry <= a_trueBdry[:,1,:,0]

        #print('here')
        #print(bdry_upperCheck_trueBdry.shape)
        #print(stat_FcBdry.shape)
        #print(a_trueBdry[:,1,:,0].shape)

        # -------------------------------------------------------------------
        # Work out whether simulation observed successful sets.
        # -------------------------------------------------------------------
        # Record if we saw a violation in the true boundary based sets
        trueBdry_success_intrp[r,:] = np.all(bdry_lowerCheck_trueBdry,axis=(1)) & np.all(bdry_upperCheck_trueBdry,axis=(1))# MARKER: AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success_intrp[r,:] = np.all(bdry_lowerCheck_estBdry,axis=(1)) & np.all(bdry_upperCheck_estBdry,axis=(1)) # MARKER: AXES WONT WORK FOR 3D ATM


    # For the interpolated boundary success checks, we still need to do the 
    # voxelwise checks as well. This will take care of that.
    trueBdry_success_intrp = trueBdry_success_intrp*trueBdry_success
    estBdry_success_intrp = estBdry_success_intrp*estBdry_success

    # Coverage probabilities
    coverage_trueBdry = np.mean(trueBdry_success,axis=0)
    coverage_estBdry = np.mean(estBdry_success,axis=0)

    # Coverage probabilities
    coverage_trueBdry_intrp = np.mean(trueBdry_success_intrp,axis=0)
    coverage_estBdry_intrp = np.mean(estBdry_success_intrp,axis=0)

    print('Coverage: ', coverage_estBdry_intrp)
    print('Coverage: ', coverage_trueBdry_intrp)

    # Make results folder
    if not os.path.exists(os.path.join(simDir, 'RawResults')):
        os.mkdir(os.path.join(simDir, 'RawResults'))

    # Save the violations to a file
    append_to_file(os.path.join(simDir, 'RawResults', 'trueSuccess'+str(nSub)+'.csv'), trueBdry_success) # Successes based on the true boundary (assessed without interpolation)
    append_to_file(os.path.join(simDir, 'RawResults', 'estSuccess'+str(nSub)+'.csv'), estBdry_success) # Successes based on the interpolated boundary (assessed without interpolation) 
    append_to_file(os.path.join(simDir, 'RawResults', 'trueSuccess'+str(nSub)+'_intrp.csv'), trueBdry_success_intrp) # Successes based on the true boundary (assessed with interpolation)
    append_to_file(os.path.join(simDir, 'RawResults', 'estSuccess'+str(nSub)+'_intrp.csv'), estBdry_success_intrp) # Successes based on the interpolated boundary (assessed with interpolation) 

    # Save the computation times
    t2overall = time.time()
    append_to_file(os.path.join(simDir, 'RawResults', 'computationTime'+str(nSub)+'.csv'), np.array([t2overall-t1overall]))

#SpatialSims_2mu('/home/tommaullin/Documents/ConfSets/sim1/cfgs/cfg623.yml')