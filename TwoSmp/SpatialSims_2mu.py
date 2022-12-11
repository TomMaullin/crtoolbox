import os
import time
import numpy as np
from lib.generateData import *
from lib.boundary import *
from lib.fileio import *
import yaml
from scipy.ndimage.measurements import label
import matplotlib.pyplot as plt

# ===========================================================================
#
# Simulation where `a' is controlled by
#
#       a = max(sup_(dAc1\intersectAc2) |g1|, sup_(Ac1\intersectAc2d)|g2|)
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
#                  shortening for intersection and union so I have used F and
#                  G for shorthand for intersection and union, respectively,
#                  wherever they occur. FsdG is used in this to represent
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

    # If we're running the method seperately for each field and intersecting
    # the results (for sims 29 and 30), then we move to the 'seperate'
    # function.
    if 'Seperate' in inputs:
        if inputs['Seperate']:
            SpatialSims_2mu_seperate(ipath)
            return(None)

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

    # Get noise spec
    noiseSpec1 = inputs['noise1']

    # Reforemat FWHM for the noise in field 1
    noiseSpec1['FWHM'] = str2vec(noiseSpec1['FWHM']) 

    # Get noise spec
    noiseSpec2 = inputs['noise2']

    # Reforemat FWHM for the noise in field 2
    noiseSpec2['FWHM'] = str2vec(noiseSpec2['FWHM']) 

    # Correlation between noise fields
    if 'noiseCorr' in inputs:
        noiseCorr = np.float(inputs['noiseCorr'])
    else:
        noiseCorr = None

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


    # Is there some signal we want both images to have?
    if 'muSpecBoth' in inputs:

        # Get specification for mu overlap
        muSpecBoth = inputs['muBoth']

        # Reformat center and fwhm if needed
        if 'center' in muSpecBoth:
            muSpecBoth['center']=eval(muSpecBoth['center'])
        if 'fwhm' in muSpecBoth:
            muSpecBoth['fwhm']=eval(muSpecBoth['fwhm'])

    else:

        # Otherwise set it to none
        muSpecBoth = None

    # ID for the configuration
    cfgId = inputs['cfgId']

    # Bootstrap mode
    mode = inputs['mode']

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
    boot_dim = np.array([nSub, 1]) # : COULD MAKE  np.array([batchBoot, nSub, 1])


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

    # Initialize for saving times
    times = np.zeros((nReals,1))

    # Loop through realizations
    for r in np.arange(nReals):

        print('r: ', r)
        # -------------------------------------------------------------------
        # Data generation
        # -------------------------------------------------------------------

        # Obtain data
        data1, data2, mu1, mu2 = get_data(muSpec1,muSpec2,noiseSpec1,noiseSpec2, data_dim, noiseCorr, muSpecBoth)

        #print('data shapes: ', data1.shape, data2.shape, mu1.shape, mu2.shape)

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
        resid1_dFc_concat = get_bdry_values_concat(resid1, Fc_bdry_locs)
        resid2_dFc_concat = get_bdry_values_concat(resid2, Fc_bdry_locs)

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


        # print('shapes here : ', mu1_Fc_bdry_concat.shape,resid1_dFc_concat.shape)

        # Get locations where outer mu1 and mu2 are greater than c
        d1Fc_loc = np.where(mu2_Fc_bdry_concat[0,:,1]>c)[0]
        d2Fc_loc = np.where(mu1_Fc_bdry_concat[0,:,1]>c)[0]
        d12Fc_loc = np.where((mu2_Fc_bdry_concat[0,:,1]<= c)*(mu1_Fc_bdry_concat[0,:,1]<=c))[0]

        # print(resid1_dFc_concat[:,d1Fc_loc,:].shape)
        # print(resid2_dFc_concat[:,d2Fc_loc,:].shape)

        # Obtain MuHat along Fc
        muHat1_FcHat_bdry_concat = get_bdry_values_concat(muHat1, FcHat_bdry_locs)
        muHat2_FcHat_bdry_concat = get_bdry_values_concat(muHat2, FcHat_bdry_locs)

        # print('shapes here : ', muHat1_FcHat_bdry_concat.shape,resid1_FcHat_bdry_concat.shape)

        # Get locations where outer muHat1 and muHat2 are greater than c
        d1FcHat_loc = np.where(muHat2_FcHat_bdry_concat[0,:,1]>c)[0]
        d2FcHat_loc = np.where(muHat1_FcHat_bdry_concat[0,:,1]>c)[0]
        d12FcHat_loc = np.where((muHat2_FcHat_bdry_concat[0,:,1]<= c)*(muHat1_FcHat_bdry_concat[0,:,1]<=c))[0]

        # print(resid1_FcHat_bdry_concat[:,d1FcHat_loc,:].shape)
        # print(resid2_FcHat_bdry_concat[:,d2FcHat_loc,:].shape)

        # -------------------------------------------------------------------
        # Residuals along dkFc and dkFcHat boundaries (Array version)
        # -------------------------------------------------------------------
        # In this simulation we are bootstrapping the residuals for field 1
        # along dAc1 intersect Ac2 i.e. along dF where mu2 > c. And vice versa
        # for field 2.
        resid1_d1Fc_concat = resid1_dFc_concat[:,d1Fc_loc,:]
        resid2_d2Fc_concat = resid2_dFc_concat[:,d2Fc_loc,:]

        # In this simulation we are bootstrapping the residuals for field 1
        # along dAcHat1 intersect AcHat2 i.e. along dF where muHat2 > c. And 
        # vice versa for field 2.
        resid1_d1FcHat_concat = resid1_FcHat_bdry_concat[:,d1FcHat_loc,:]
        resid2_d2FcHat_concat = resid2_FcHat_bdry_concat[:,d2FcHat_loc,:]

        # If we are in mode 2 or 3 we need to bootstrap along the intersection
        # boundary as well
        if mode == 2 or mode == 3:

            # For true boundary
            resid1_d12Fc_concat = resid1_dFc_concat[:,d12Fc_loc,:]
            resid2_d12Fc_concat = resid2_dFc_concat[:,d12Fc_loc,:]

            # For estimated boundary
            resid1_d12FcHat_concat = resid1_FcHat_bdry_concat[:,d12FcHat_loc,:]
            resid2_d12FcHat_concat = resid2_FcHat_bdry_concat[:,d12FcHat_loc,:]

        # -------------------------------------------------------------------
        # Mu and MuHat along dkFc and dkFcHat boundaries (Array version)
        # -------------------------------------------------------------------
        # In this simulation we are bootstrapping the residuals for field 1
        # along dAc1 intersect Ac2 i.e. along dF where mu2 > c. And vice versa
        # for field 2.
        mu1_d1Fc_concat = mu1_Fc_bdry_concat[:,d1Fc_loc,:]
        mu2_d2Fc_concat = mu2_Fc_bdry_concat[:,d2Fc_loc,:]

        # In this simulation we are bootstrapping the residuals for field 1
        # along dAcHat1 intersect AcHat2 i.e. along dF where muHat2 > c. And 
        # vice versa for field 2.
        muHat1_d1FcHat_concat = muHat1_FcHat_bdry_concat[:,d1FcHat_loc,:]
        muHat2_d2FcHat_concat = muHat2_FcHat_bdry_concat[:,d2FcHat_loc,:]

        # If we are in mode 2 or 3 we need to bootstrap along the intersection
        # boundary as well
        if mode == 2 or mode == 3:
            
            # For true boundary
            mu1_d12Fc_concat = mu1_Fc_bdry_concat[:,d12Fc_loc,:]
            mu2_d12Fc_concat = mu2_Fc_bdry_concat[:,d12Fc_loc,:]

            # For estimated boundary
            muHat1_d12FcHat_concat = muHat1_FcHat_bdry_concat[:,d12FcHat_loc,:]
            muHat2_d12FcHat_concat = muHat2_FcHat_bdry_concat[:,d12FcHat_loc,:]


        # print('marker')
        # print(mu1_d1Fc_concat.shape)
        # print(mu1_Fc_bdry_concat.shape)
        # print(d2Fc_loc.shape)
        # print(mu2_d2Fc_concat.shape)
        # print(mu2_Fc_bdry_concat.shape)
        # print(d1Fc_loc.shape)
        # print(muHat1_d1FcHat_concat.shape)

        # print('marker2')
        # print(d1FcHat_loc.shape)
        # print(d2FcHat_loc.shape)
        # print(d12FcHat_loc.shape)
        # print(muHat1_FcHat_bdry_concat.shape)

        # -------------------------------------------------------------------
        # Interpolation weights dkFc and dkFcHat boundaries (Array version)
        # -------------------------------------------------------------------
        # Obtain the weights along the boundary for dkFc
        d1Fc_bdry_weights_concat = get_bdry_weights_concat(mu1_d1Fc_concat, c)
        d2Fc_bdry_weights_concat = get_bdry_weights_concat(mu2_d2Fc_concat, c)

        # Obtain the weights along the boundary for dkFcHat
        d1FcHat_weights_concat = get_bdry_weights_concat(muHat1_d1FcHat_concat, c)
        d2FcHat_weights_concat = get_bdry_weights_concat(muHat2_d2FcHat_concat, c)

        # If we are in mode 2 or 3 we need to bootstrap along the intersection
        # boundary as well
        if mode == 2 or mode == 3:
            
            # For true boundary
            d12Fc_mu1_bdry_weights_concat = get_bdry_weights_concat(mu1_d12Fc_concat, c)
            d12Fc_mu2_bdry_weights_concat = get_bdry_weights_concat(mu2_d12Fc_concat, c)

            # For estimated boundary
            d12FcHat_muHat1_weights_concat = get_bdry_weights_concat(muHat1_d12FcHat_concat, c)
            d12FcHat_muHat2_weights_concat = get_bdry_weights_concat(muHat2_d12FcHat_concat, c)

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
        max_g1_d1Fc = np.zeros(nBoot)
        max_g2_d2Fc = np.zeros(nBoot)
        max_g1_d1FcHat = np.zeros(nBoot)
        max_g2_d2FcHat = np.zeros(nBoot)

        if mode == 2:

            # Empty bootstrap store for intersection boundary
            max_g1g2_d12Fc = np.zeros(nBoot)
            max_g1g2_d12FcHat = np.zeros(nBoot)

        if mode == 3:

            # Empty bootstrap store for intersection boundary
            min_g1g2_d12Fc = np.zeros(nBoot)
            min_g1g2_d12FcHat = np.zeros(nBoot)

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

            # Bootstrap residuals along d1Fc and d2Fc
            boot_resid1_d1Fc_concat = boot_vars*resid1_d1Fc_concat
            boot_resid2_d2Fc_concat = boot_vars*resid2_d2Fc_concat

            # Bootstrap residuals along d1FcHat and d2FcHat
            boot_resid1_d1FcHat_concat = boot_vars*resid1_d1FcHat_concat
            boot_resid2_d2FcHat_concat = boot_vars*resid2_d2FcHat_concat

            # If we are in mode 2 or 3 we need to bootstrap along the intersection
            # boundary as well
            if mode == 2 or mode == 3:
            
                # Bootstrap residuals along d12Fc
                boot_resid1_d12Fc_concat = boot_vars*resid1_d12Fc_concat
                boot_resid2_d12Fc_concat = boot_vars*resid2_d12Fc_concat

                # Bootstrap residuals along d12FcHat
                boot_resid1_d12FcHat_concat = boot_vars*resid1_d12FcHat_concat
                boot_resid2_d12FcHat_concat = boot_vars*resid2_d12FcHat_concat

            # ------------------------------------------------------------------------------------------
            # True boundary
            # ------------------------------------------------------------------------------------------

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of d1Fc. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g1_d1Fc_concat = np.zeros(boot_resid1_d1Fc_concat.shape[-2:])
            boot_g1_d1Fc_concat[...,0] = np.sum(boot_resid1_d1Fc_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g1_d1Fc_concat[...,1] = np.sum(boot_resid1_d1Fc_concat[...,1], axis=0)/np.sqrt(nSub)

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of d2Fc. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g2_d2Fc_concat = np.zeros(boot_resid2_d2Fc_concat.shape[-2:])
            boot_g2_d2Fc_concat[...,0] = np.sum(boot_resid2_d2Fc_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g2_d2Fc_concat[...,1] = np.sum(boot_resid2_d2Fc_concat[...,1], axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along d1Fc. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma1_boot_d1Fc_concat = np.zeros(boot_resid1_d1Fc_concat.shape[-2:])
            sigma1_boot_d1Fc_concat[...,0] = np.std(boot_resid1_d1Fc_concat[...,0], axis=0, ddof=1)
            sigma1_boot_d1Fc_concat[...,1] = np.std(boot_resid1_d1Fc_concat[...,1], axis=0, ddof=1)

            # Obtain bootstrap standard deviations along d2Fc. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma2_boot_d2Fc_concat = np.zeros(boot_resid2_d2Fc_concat.shape[-2:])
            sigma2_boot_d2Fc_concat[...,0] = np.std(boot_resid2_d2Fc_concat[...,0], axis=0, ddof=1)
            sigma2_boot_d2Fc_concat[...,1] = np.std(boot_resid2_d2Fc_concat[...,1], axis=0, ddof=1)

            # Divide by the boostrap standard deviation on d1Fc
            boot_g1_d1Fc_concat = boot_g1_d1Fc_concat/sigma1_boot_d1Fc_concat

            # Divide by the boostrap standard deviation on d2Fc
            boot_g2_d2Fc_concat = boot_g2_d2Fc_concat/sigma2_boot_d2Fc_concat

            # ------------------------------------------------------------------------------------------
            # True boundary (mode == 2 or 3)
            # ------------------------------------------------------------------------------------------

            # If we are in mode 2 or 3 we need to bootstrap along the intersection
            # boundary as well
            if mode == 2 or mode == 3:

                # Sum across subjects to get the bootstrapped a values along
                # the boundary of d12Fc. (Note: For some reason this is 
                # much faster if performed seperately for each of the last rows. 
                # I am still looking into why this is)
                boot_g1_d12Fc_concat = np.zeros(boot_resid1_d12Fc_concat.shape[-2:])
                boot_g1_d12Fc_concat[...,0] = np.sum(boot_resid1_d12Fc_concat[...,0], axis=0)/np.sqrt(nSub)
                boot_g1_d12Fc_concat[...,1] = np.sum(boot_resid1_d12Fc_concat[...,1], axis=0)/np.sqrt(nSub)

                # Sum across subjects to get the bootstrapped a values along
                # the boundary of d12Fc. (Note: For some reason this is 
                # much faster if performed seperately for each of the last rows. 
                # I am still looking into why this is)
                boot_g2_d12Fc_concat = np.zeros(boot_resid2_d12Fc_concat.shape[-2:])
                boot_g2_d12Fc_concat[...,0] = np.sum(boot_resid2_d12Fc_concat[...,0], axis=0)/np.sqrt(nSub)
                boot_g2_d12Fc_concat[...,1] = np.sum(boot_resid2_d12Fc_concat[...,1], axis=0)/np.sqrt(nSub)

                # Obtain bootstrap standard deviations along d12Fc. (Note: For some reason this is 
                # much faster if performed seperately for each of the last rows. I am still looking
                # into why this is)
                sigma1_boot_d12Fc_concat = np.zeros(boot_resid1_d12Fc_concat.shape[-2:])
                sigma1_boot_d12Fc_concat[...,0] = np.std(boot_resid1_d12Fc_concat[...,0], axis=0, ddof=1)
                sigma1_boot_d12Fc_concat[...,1] = np.std(boot_resid1_d12Fc_concat[...,1], axis=0, ddof=1)

                # Obtain bootstrap standard deviations along d12Fc. (Note: For some reason this is 
                # much faster if performed seperately for each of the last rows. I am still looking
                # into why this is)
                sigma2_boot_d12Fc_concat = np.zeros(boot_resid2_d12Fc_concat.shape[-2:])
                sigma2_boot_d12Fc_concat[...,0] = np.std(boot_resid2_d12Fc_concat[...,0], axis=0, ddof=1)
                sigma2_boot_d12Fc_concat[...,1] = np.std(boot_resid2_d12Fc_concat[...,1], axis=0, ddof=1)

                # Divide by the boostrap standard deviation on d12Fc
                boot_g1_d12Fc_concat = boot_g1_d12Fc_concat/sigma1_boot_d12Fc_concat

                # Divide by the boostrap standard deviation on d12Fc
                boot_g2_d12Fc_concat = boot_g2_d12Fc_concat/sigma2_boot_d12Fc_concat

            # ------------------------------------------------------------------------------------------
            # Estimated boundary
            # ------------------------------------------------------------------------------------------

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of d1FcHat. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g1_d1FcHat_concat = np.zeros(boot_resid1_d1FcHat_concat.shape[-2:])
            boot_g1_d1FcHat_concat[...,0] = np.sum(boot_resid1_d1FcHat_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g1_d1FcHat_concat[...,1] = np.sum(boot_resid1_d1FcHat_concat[...,1], axis=0)/np.sqrt(nSub)

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of d2FcHat. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g2_d2FcHat_concat = np.zeros(boot_resid2_d2FcHat_concat.shape[-2:])
            boot_g2_d2FcHat_concat[...,0] = np.sum(boot_resid2_d2FcHat_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g2_d2FcHat_concat[...,1] = np.sum(boot_resid2_d2FcHat_concat[...,1], axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along d1FcHat. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma1_boot_d1FcHat_concat = np.zeros(boot_resid1_d1FcHat_concat.shape[-2:])
            sigma1_boot_d1FcHat_concat[...,0] = np.std(boot_resid1_d1FcHat_concat[...,0], axis=0, ddof=1)
            sigma1_boot_d1FcHat_concat[...,1] = np.std(boot_resid1_d1FcHat_concat[...,1], axis=0, ddof=1)

            # Obtain bootstrap standard deviations along d2FcHat. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma2_boot_d2FcHat_concat = np.zeros(boot_resid2_d2FcHat_concat.shape[-2:])
            sigma2_boot_d2FcHat_concat[...,0] = np.std(boot_resid2_d2FcHat_concat[...,0], axis=0, ddof=1)
            sigma2_boot_d2FcHat_concat[...,1] = np.std(boot_resid2_d2FcHat_concat[...,1], axis=0, ddof=1)

            # Divide by the boostrap standard deviation on d1FcHat
            boot_g1_d1FcHat_concat = boot_g1_d1FcHat_concat/sigma1_boot_d1FcHat_concat

            # Divide by the boostrap standard deviation on d2FcHat
            boot_g2_d2FcHat_concat = boot_g2_d2FcHat_concat/sigma2_boot_d2FcHat_concat
                
            # ------------------------------------------------------------------------------------------
            # Estimated boundary (mode == 2 or 3)
            # ------------------------------------------------------------------------------------------

            # If we are in mode 2 or 3 we need to bootstrap along the intersection
            # boundary as well
            if mode == 2 or mode == 3:

                # Sum across subjects to get the bootstrapped a values along
                # the boundary of d12FcHat. (Note: For some reason this is 
                # much faster if performed seperately for each of the last rows. 
                # I am still looking into why this is)
                boot_g1_d12FcHat_concat = np.zeros(boot_resid1_d12FcHat_concat.shape[-2:])
                boot_g1_d12FcHat_concat[...,0] = np.sum(boot_resid1_d12FcHat_concat[...,0], axis=0)/np.sqrt(nSub)
                boot_g1_d12FcHat_concat[...,1] = np.sum(boot_resid1_d12FcHat_concat[...,1], axis=0)/np.sqrt(nSub)

                # Sum across subjects to get the bootstrapped a values along
                # the boundary of d12FcHat. (Note: For some reason this is 
                # much faster if performed seperately for each of the last rows. 
                # I am still looking into why this is)
                boot_g2_d12FcHat_concat = np.zeros(boot_resid2_d12FcHat_concat.shape[-2:])
                boot_g2_d12FcHat_concat[...,0] = np.sum(boot_resid2_d12FcHat_concat[...,0], axis=0)/np.sqrt(nSub)
                boot_g2_d12FcHat_concat[...,1] = np.sum(boot_resid2_d12FcHat_concat[...,1], axis=0)/np.sqrt(nSub)

                # Obtain bootstrap standard deviations along d12FcHat. (Note: For some reason this is 
                # much faster if performed seperately for each of the last rows. I am still looking
                # into why this is)
                sigma1_boot_d12FcHat_concat = np.zeros(boot_resid1_d12FcHat_concat.shape[-2:])
                sigma1_boot_d12FcHat_concat[...,0] = np.std(boot_resid1_d12FcHat_concat[...,0], axis=0, ddof=1)
                sigma1_boot_d12FcHat_concat[...,1] = np.std(boot_resid1_d12FcHat_concat[...,1], axis=0, ddof=1)

                # Obtain bootstrap standard deviations along d12FcHat. (Note: For some reason this is 
                # much faster if performed seperately for each of the last rows. I am still looking
                # into why this is)
                sigma2_boot_d12FcHat_concat = np.zeros(boot_resid2_d12FcHat_concat.shape[-2:])
                sigma2_boot_d12FcHat_concat[...,0] = np.std(boot_resid2_d12FcHat_concat[...,0], axis=0, ddof=1)
                sigma2_boot_d12FcHat_concat[...,1] = np.std(boot_resid2_d12FcHat_concat[...,1], axis=0, ddof=1)

                # Divide by the boostrap standard deviation on d12FcHat
                boot_g1_d12FcHat_concat = boot_g1_d12FcHat_concat/sigma1_boot_d12FcHat_concat

                # Divide by the boostrap standard deviation on d12FcHat
                boot_g2_d12FcHat_concat = boot_g2_d12FcHat_concat/sigma2_boot_d12FcHat_concat

            # ------------------------------------------------------------------------------------------
            # Interpolation
            # ------------------------------------------------------------------------------------------

            # Interpolation for d1Fc and d2Fc boundary
            boot_g1_d1Fc_concat = get_bdry_vals_interpolated_concat(boot_g1_d1Fc_concat,d1Fc_bdry_weights_concat)
            boot_g2_d2Fc_concat = get_bdry_vals_interpolated_concat(boot_g2_d2Fc_concat,d2Fc_bdry_weights_concat)

            # Interpolation for d1FcHat and d2FcHat boundary
            boot_g1_d1FcHat_concat = get_bdry_vals_interpolated_concat(boot_g1_d1FcHat_concat,d1FcHat_weights_concat)
            boot_g2_d2FcHat_concat = get_bdry_vals_interpolated_concat(boot_g2_d2FcHat_concat,d2FcHat_weights_concat)

            # If we are in mode 2 or 3 we need to bootstrap along the intersection
            # boundary as well
            if mode == 2 or mode == 3:

                # Interpolation for d1Fc and d2Fc boundary
                boot_g1_d12Fc_concat = get_bdry_vals_interpolated_concat(boot_g1_d12Fc_concat,d12Fc_mu1_bdry_weights_concat)
                boot_g2_d12Fc_concat = get_bdry_vals_interpolated_concat(boot_g2_d12Fc_concat,d12Fc_mu2_bdry_weights_concat)

                # Interpolation for d12FcHat and d12FcHat boundary
                boot_g1_d12FcHat_concat = get_bdry_vals_interpolated_concat(boot_g1_d12FcHat_concat,d12FcHat_muHat1_weights_concat)
                boot_g2_d12FcHat_concat = get_bdry_vals_interpolated_concat(boot_g2_d12FcHat_concat,d12FcHat_muHat2_weights_concat)


            # Get the maximum along d1Fc if it exists
            if np.prod(boot_g1_d1Fc_concat.shape) > 0:
                max_g1_d1Fc[b] = np.max(np.abs(boot_g1_d1Fc_concat)) 


            # Get the maximum along d2Fc if it exists
            if np.prod(boot_g2_d2Fc_concat.shape) > 0:
                max_g2_d2Fc[b] = np.max(np.abs(boot_g2_d2Fc_concat)) 

            # Get the maximum along d1FcHat if it exists
            if np.prod(boot_g1_d1FcHat_concat.shape) > 0:
                max_g1_d1FcHat[b] = np.max(np.abs(boot_g1_d1FcHat_concat)) 


            # Get the maximum along d2FcHat if it exists
            if np.prod(boot_g2_d2FcHat_concat.shape) > 0:
                max_g2_d2FcHat[b] = np.max(np.abs(boot_g2_d2FcHat_concat)) 

            if mode == 2:

                # Get the maximum along d12Fc if it exists
                if np.prod(boot_g1_d12Fc_concat.shape) > 0:

                    # Maximum of both |g1| and |g2| on d12Fc boundary
                    max_g1g2_d12Fc[b] = np.maximum(np.max(np.abs(boot_g1_d12Fc_concat)),np.max(np.abs(boot_g2_d12Fc_concat)))

                # Get the maximum along d12FcHat if it exists
                if np.prod(boot_g1_d12FcHat_concat.shape) > 0:

                    # Maximum of both |g1| and |g2| on d12Fc boundary
                    max_g1g2_d12FcHat[b] = np.maximum(np.max(np.abs(boot_g1_d12FcHat_concat)),np.max(np.abs(boot_g2_d12FcHat_concat)))

            if mode == 3:

                #print('marker4')

                # Get the maximum along d12Fc if it exists
                if np.prod(boot_g1_d12Fc_concat.shape) > 0:

                    #print('marker3')

                    # Minimum of both g1 and g2 on d12Fc boundary (note: the absolute values must be outside
                    # the minimum in this mode, unlike in mode 2 where they were inside)
                    min_g1g2_d12Fc[b] = np.max(np.abs(np.minimum(boot_g1_d12Fc_concat,boot_g2_d12Fc_concat)))

                    # print('1: ', min_g1g2_d12Fc)
                    # print(min_g1g2_d12Fc[b])
 
                # Get the maximum along d12FcHat if it exists
                if np.prod(boot_g1_d12FcHat_concat.shape) > 0:
                    
                    # Minimum of both g1 and g2 on d12Fc boundary (note: the absolute values must be outside
                    # the minimum in this mode, unlike in mode 2 where they were inside)
                    min_g1g2_d12FcHat[b] = np.max(np.abs(np.minimum(boot_g1_d12FcHat_concat,boot_g2_d12FcHat_concat)))


            # -------------------------------------------------------------------------------------


        # In mode 1, we only bootstrap G1 along d1Fc and G2 along d2Fc
        if mode == 1:

            # Get the maximum needed for the true boundary
            max_g_dFc = np.maximum(max_g1_d1Fc,max_g2_d2Fc)

            # Get the maximum needed for the estimated boundary
            max_g_dFcHat = np.maximum(max_g1_d1FcHat,max_g2_d2FcHat)

        # In mode 2, we bootstrap G1 along d1Fc union d12Fc and G2 along d2Fc union d12Fc
        if mode == 2:

            # Get the maximum needed for the true boundary
            max_g_dFc = np.maximum(np.maximum(max_g1_d1Fc,max_g2_d2Fc),max_g1g2_d12Fc)

            # Get the maximum needed for the estimated boundary
            max_g_dFcHat = np.maximum(np.maximum(max_g1_d1FcHat,max_g2_d2FcHat),max_g1g2_d12FcHat)


        # In mode 3, we bootstrap G1 along d1Fc and G2 along d2Fc and min(G1,G2) along d12Fc
        if mode == 3:

            # print('3: ', min_g1g2_d12Fc)
            # Get the maximum needed for the true boundary
            max_g_dFc = np.maximum(np.maximum(max_g1_d1Fc,max_g2_d2Fc),min_g1g2_d12Fc)

            # Get the maximum needed for the estimated boundary
            max_g_dFcHat = np.maximum(np.maximum(max_g1_d1FcHat,max_g2_d2FcHat),min_g1g2_d12FcHat)

            # print('active')
            # print(max_g1_d1Fc.shape,max_g2_d2Fc.shape,min_g1g2_d12Fc.shape)
            # print('1')
            # print(max_g1_d1Fc)
            # print('2')
            # print(max_g2_d2Fc)
            # print('1and2')
            # print('4: ', min_g1g2_d12Fc)
            # print('max')
            # print(max_g_dFc)

            # print('active2')
            # print(max_g1_d1FcHat.shape,max_g2_d2FcHat.shape,min_g1g2_d12FcHat.shape)
            # print(max_g_dFcHat)
            
        t2 = time.time()
        # print('Bootstrap time: ', t2-t1)

        # -------------------------------------------------------------------
        # Obtaining a from percentiles of the max distribution
        # -------------------------------------------------------------------

        # Drop the instances where the boundary length was zero
        max_g_dFc = max_g_dFc[max_g_dFc!=0]
        max_g_dFcHat = max_g_dFcHat[max_g_dFcHat!=0]

        # If we have recorded values get their quantiles
        if (np.prod(max_g_dFc.shape) > 0):
            # Get the a estimates for the true boundary
            a_trueBdry = np.percentile(max_g_dFc, 100*p).reshape(nPvals,1,1,1)
        else:
            # Set to inf by default
            a_trueBdry = np.Inf*np.ones((nPvals,1,1,1))


        # If we have recorded values get their quantiles
        if (np.prod(max_g_dFcHat.shape) > 0):
            # Get the a estimates for the estimated boundary
            a_estBdry = np.percentile(max_g_dFcHat, 100*p).reshape(nPvals,1,1,1) # [pvals, 1, [1 for _ in dim>1]]
        else:
            # Set to inf by default
            a_estBdry = np.Inf*np.ones((nPvals,1,1,1))

        # Reformat them to an array form useful for boolean operation
        a_trueBdry = np.concatenate((-a_trueBdry,a_trueBdry),axis=1)
        a_estBdry = np.concatenate((-a_estBdry,a_estBdry),axis=1)

        # print('a')
        # print(a_trueBdry)
        # print(a_estBdry)

        # Get the statistic field which defined Achat^{+/-,1/2}
        g1 = ((muHat1-c)/(sigma1*tau))
        g2 = ((muHat2-c)/(sigma2*tau))

        # Minimum for intersection
        stat = np.minimum(g1.reshape(1,(*muHat1.shape)),g2.reshape(1,(*muHat1.shape)))
        stat = stat.reshape(stat.shape[-2],stat.shape[-1])

        # Obtain FcHat^+ and FcHat^- based on a from the true boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        FcHat_pm_trueBdry = stat >= a_trueBdry

        # Obtain FcHat^+ and FcHat^- based on a from the estimated boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        FcHat_pm_estBdry = stat >= a_estBdry


        # End timer (we have the confidence sets now)
        t2 = time.time()

        # Save time 
        times[r,:] = t2-t1
        
        # -------------------------------------------------------------------
        # Images
        # -------------------------------------------------------------------
        if r==1 and inputs['figGen']:

            # Indices for 0.8,0.9 and 0.95
            pInds = [16,18,19]

            for pInd in pInds:

                # P value of interest
                pVal = p[pInd]

                # FcHat plus (based on estimated bdry)
                plt.figure(16)
                plt.imshow(1*FcHat_pm_estBdry[pInd,0,...])
                plt.savefig(os.path.join(figDir, 'FcHat_plus_estBdry_p'+str(int(100*pVal))+'_cfg'+str(cfgId)+'.png'))

                # FcHat minus (based on estimated bdry)
                plt.figure(17)
                plt.imshow(1*FcHat_pm_estBdry[pInd,1,...])
                plt.savefig(os.path.join(figDir, 'FcHat_minus_estBdry_p'+str(int(100*pVal))+'_cfg'+str(cfgId)+'.png'))

                # FcHat plus (based on true bdry)
                plt.figure(18)
                plt.imshow(1*FcHat_pm_trueBdry[pInd,0,...])
                plt.savefig(os.path.join(figDir, 'FcHat_plus_trueBdry_p'+str(int(100*pVal))+'_cfg'+str(cfgId)+'.png'))

                # FcHat minus (based on true bdry)
                plt.figure(19)
                plt.imshow(1*FcHat_pm_trueBdry[pInd,1,...])
                plt.savefig(os.path.join(figDir, 'FcHat_minus_trueBdry_p'+str(int(100*pVal))+'_cfg'+str(cfgId)+'.png'))


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

        # -------------------------------------------------------------------
        # Check whether there were any boundary violations using voxelwise
        # set logic (checking if voxels existed in one set but not another, 
        # etc)
        # -------------------------------------------------------------------
        #print(a_trueBdry)
        # # FcHat plus (based on estimated bdry)
        # plt.figure(0)
        # plt.imshow(1*FcHat_pm_trueBdry[0,1,...])
        # plt.colorbar()

        # # FcHat minus (based on estimated bdry)
        # plt.figure(1)
        # plt.imshow(1*(~Fc[0,...]))
        # plt.colorbar()
        # plt.show()
        # plt.clf()


        #print('anys: ', np.any(FcHatp_sub_Fc_trueBdry,axis=(1,2)), np.any(Fc_sub_FcHatm_trueBdry,axis=(1,2)))
        # Record if we saw a violation in the true boundary based sets
        trueBdry_success[r,:] = 1-(np.any(FcHatp_sub_Fc_trueBdry,axis=(1,2)) | np.any(Fc_sub_FcHatm_trueBdry,axis=(1,2))) # : AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success[r,:] = 1-(np.any(FcHatp_sub_Fc_estBdry,axis=(1,2)) | np.any(Fc_sub_FcHatm_estBdry,axis=(1,2)))# : AXES WONT WORK FOR 3D ATM

        # -------------------------------------------------------------------
        # Get stat along the Fc boundary
        # -------------------------------------------------------------------


        # ==================================================================================================================

        # Get locations where outer mu1 and mu2 are greater than c
        d1Fc_loc = np.where((mu2_Fc_bdry_concat[0,:,1]>c))[0]
        d2Fc_loc = np.where((mu1_Fc_bdry_concat[0,:,1]>c))[0]
        d12Fc_loc = np.where((mu2_Fc_bdry_concat[0,:,1]<= c)*(mu1_Fc_bdry_concat[0,:,1]<=c))[0]



        # Obtain g1 and g2 along the boundary for Fc
        g1_dFc_concat = get_bdry_values_concat(g1, Fc_bdry_locs)
        g2_dFc_concat = get_bdry_values_concat(g2, Fc_bdry_locs)


        # g1 and g2 along d1Fc
        g1_d1Fc_concat = g1_dFc_concat[:,d1Fc_loc,:]
        g2_d1Fc_concat = g2_dFc_concat[:,d1Fc_loc,:]

        # g1 and g2 along d2Fc
        g1_d2Fc_concat = g1_dFc_concat[:,d2Fc_loc,:]
        g2_d2Fc_concat = g2_dFc_concat[:,d2Fc_loc,:]

        # g1 and g2 along d12Fc
        g1_d12Fc_concat = g1_dFc_concat[:,d12Fc_loc,:]
        g2_d12Fc_concat = g2_dFc_concat[:,d12Fc_loc,:]

        # Interpolation for d1Fc boundary (we perform this interpolation based on the weights of
        # mu1 as we are on the boundary of Ac1)
        g1_d1Fc_concat = get_bdry_vals_interpolated_concat(g1_d1Fc_concat,d1Fc_bdry_weights_concat)
        g2_d1Fc_concat = get_bdry_vals_interpolated_concat(g2_d1Fc_concat,d1Fc_bdry_weights_concat)


        # Interpolation for d2Fc boundary (we perform this interpolation based on the weights of
        # mu2 as we are on the boundary of Ac2)
        g1_d2Fc_concat = get_bdry_vals_interpolated_concat(g1_d2Fc_concat,d2Fc_bdry_weights_concat)
        g2_d2Fc_concat = get_bdry_vals_interpolated_concat(g2_d2Fc_concat,d2Fc_bdry_weights_concat)


        # If we are in mode 2 or 3 we need to interpolate along the intersection
        # boundary as well
        if mode == 2 or mode == 3:

            # Interpolation for d1Fc and d2Fc boundary
            g1_d12Fc_concat = get_bdry_vals_interpolated_concat(g1_d12Fc_concat,d12Fc_mu1_bdry_weights_concat)
            g2_d12Fc_concat = get_bdry_vals_interpolated_concat(g2_d12Fc_concat,d12Fc_mu2_bdry_weights_concat)


        # Get minimum along interpolated d1Fc boundary
        ming1g2_d1Fc_concat = np.minimum(g1_d1Fc_concat,g2_d1Fc_concat)

        # Get minimum along interpolated d2Fc boundary
        ming1g2_d2Fc_concat = np.minimum(g1_d2Fc_concat,g2_d2Fc_concat)

        # Get minimum along interpolated d12Fc boundary
        ming1g2_d12Fc_concat = np.minimum(g1_d12Fc_concat,g2_d12Fc_concat)

        # print('yaaah')
        # print(ming1g2_d12Fc_concat.shape)
        # print(ming1g2_d1Fc_concat)
        # print('break')
        # print(ming1g2_d2Fc_concat)
        # print(np.concatenate((ming1g2_d1Fc_concat,ming1g2_d2Fc_concat,ming1g2_d12Fc_concat),axis=-1).shape)

        # ==================================================================================================================

        # -------------------------------------------------------------------
        # Get stat along the Fc boundary
        # -------------------------------------------------------------------

        # Take the minimum of the two statistics # NTS MODE 1 OR 2 WILL CURRENTLY BREAK HERE
        stat_FcBdry = np.concatenate((ming1g2_d1Fc_concat,ming1g2_d2Fc_concat,ming1g2_d12Fc_concat),axis=-1)#stat_FcBdry.reshape(stat_FcBdry.shape[-2],stat_FcBdry.shape[-1])
        stat_FcBdry = stat_FcBdry.reshape(stat_FcBdry.shape[-2],stat_FcBdry.shape[-1])


        # # ----------------------------------------------------------------------------------------
        # # Get the values along the outer and inner boundaries
        # stat_FcBdry1 = get_bdry_values_concat(stat, Fc_bdry_locs)
        # print('MARKER 1: ', stat_FcBdry1)

        # # Interpolate to get the values along the true boundary
        # stat_FcBdry1 = get_bdry_vals_interpolated_concat(stat_FcBdry1, Fc_bdry_weights_concat)

        # # ----------------------------------------------------------------------------------------

        


        # # ----------------------------------------------------------------------------------------
        # # Get the values along the outer and inner boundaries
        # stat_FcBdry1_tmp = get_bdry_values(stat, Fc_bdry_locs)
        # print('MARKER 2: ', stat_FcBdry1)


        # # Interpolate to get the values along the true boundary
        # stat_FcBdry1_tmp = get_bdry_vals_interpolated(stat_FcBdry1_tmp, Fc_bdry_weights)
        # # ----------------------------------------------------------------------------------------

        # print(np.all(np.sort(stat_FcBdry1,axis=1)==np.sort(stat_FcBdry1_tmp,axis=1)))
        # print(stat_FcBdry1.shape)
        # print(stat_FcBdry1_tmp.shape)

        # plt.figure(0)
        # plt.hist(stat_FcBdry1.reshape(np.prod(stat_FcBdry1.shape)))#stat_FcBdry1.reshape(np.prod(stat_FcBdry.shape)))

        # plt.figure(1)
        # plt.hist(stat_FcBdry1_tmp.reshape(np.prod(stat_FcBdry1_tmp.shape)))#stat_FcBdry2.reshape(np.prod(stat_FcBdry2.shape)))
        # plt.show()

        # ming1g2_d1Fc_concat1 = stat_FcBdry1_tmp[:,d1Fc_loc]
        # ming1g2_d2Fc_concat1 = stat_FcBdry1_tmp[:,d2Fc_loc]
        # ming1g2_d12Fc_concat1 = stat_FcBdry1_tmp[:,d12Fc_loc]

        # # Take the minimum of the two statistics
        # stat_FcBdry3 = np.concatenate((ming1g2_d1Fc_concat1,ming1g2_d2Fc_concat1,ming1g2_d12Fc_concat1),axis=-1)#stat_FcBdry.reshape(stat_FcBdry.shape[-2],stat_FcBdry.shape[-1])
        # stat_FcBdry3 = stat_FcBdry3.reshape(stat_FcBdry3.shape[-2],stat_FcBdry3.shape[-1])

        # plt.figure(0)
        # plt.hist(stat_FcBdry1.reshape(np.prod(stat_FcBdry1.shape)))

        # plt.figure(1)
        # plt.hist(stat_FcBdry1_tmp.reshape(np.prod(stat_FcBdry1_tmp.shape)))

        # plt.figure(2)
        # plt.hist(ming1g2_d2Fc_concat.reshape(np.prod(ming1g2_d2Fc_concat.shape)))#stat_FcBdry1.reshape(np.prod(stat_FcBdry.shape)))

        # plt.figure(3)
        # plt.hist(ming1g2_d2Fc_concat1.reshape(np.prod(ming1g2_d2Fc_concat1.shape)))#stat_FcBdry2.reshape(np.prod(stat_FcBdry2.shape)))

        # plt.figure(4)
        # plt.hist(ming1g2_d12Fc_concat.reshape(np.prod(ming1g2_d12Fc_concat.shape)))#stat_FcBdry1.reshape(np.prod(stat_FcBdry.shape)))

        # plt.figure(5)
        # plt.hist(ming1g2_d12Fc_concat1.reshape(np.prod(ming1g2_d12Fc_concat1.shape)))#stat_FcBdry2.reshape(np.prod(stat_FcBdry2.shape)))
        # plt.show()
        # print('break')
        # print(stat_FcBdry1.reshape(np.prod(stat_FcBdry1.shape)))
        # print('break')
        # print(stat_FcBdry2.reshape(np.prod(stat_FcBdry2.shape)))
        # print('break')
        # print(a_estBdry)

        #stat_FcBdry = stat_FcBdry2

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


        #print(np.all(bdry_lowerCheck_estBdry,axis=(1)), np.all(bdry_upperCheck_estBdry,axis=(1)))

        # -------------------------------------------------------------------
        # Work out whether simulation observed successful sets.
        # -------------------------------------------------------------------
        # Record if we saw a violation in the true boundary based sets
        trueBdry_success_intrp[r,:] = np.all(bdry_lowerCheck_trueBdry,axis=(1)) & np.all(bdry_upperCheck_trueBdry,axis=(1))# : AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success_intrp[r,:] = np.all(bdry_lowerCheck_estBdry,axis=(1)) & np.all(bdry_upperCheck_estBdry,axis=(1)) # : AXES WONT WORK FOR 3D ATM


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

    # print(coverage_trueBdry)
    # print(coverage_estBdry)
    # print(coverage_trueBdry_intrp)
    # print(coverage_estBdry_intrp)

    # Make results folder
    if not os.path.exists(os.path.join(simDir, 'RawResults')):
        os.mkdir(os.path.join(simDir, 'RawResults'))

    # Save the violations to a file
    append_to_file(os.path.join(simDir, 'RawResults', 'trueSuccess.csv'), trueBdry_success) # Successes based on the true boundary (assessed without interpolation)
    append_to_file(os.path.join(simDir, 'RawResults', 'estSuccess.csv'), estBdry_success) # Successes based on the interpolated boundary (assessed without interpolation) 
    append_to_file(os.path.join(simDir, 'RawResults', 'trueSuccess_intrp.csv'), trueBdry_success_intrp) # Successes based on the true boundary (assessed with interpolation)
    append_to_file(os.path.join(simDir, 'RawResults', 'estSuccess_intrp.csv'), estBdry_success_intrp) # Successes based on the interpolated boundary (assessed with interpolation) 
    append_to_file(os.path.join(simDir, 'RawResults', 'times.csv'), times) # Times for bootstrap

    # Save the computation times
    t2overall = time.time()
    append_to_file(os.path.join(simDir, 'RawResults', 'computationTime.csv'), np.array([t2overall-t1overall]))

#SpatialSims_2mu('/home/tommaullin/Documents/ConfRes/tmp/sim15/sim15/cfgs/cfg578.yml')




# Reviewer requested - intersection of separate CRs
def SpatialSims_2mu_seperate(ipath):

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

    # Get noise spec
    noiseSpec1 = inputs['noise1']

    # Reforemat FWHM for the noise in field 1
    noiseSpec1['FWHM'] = str2vec(noiseSpec1['FWHM']) 

    # Get noise spec
    noiseSpec2 = inputs['noise2']

    # Reforemat FWHM for the noise in field 2
    noiseSpec2['FWHM'] = str2vec(noiseSpec2['FWHM']) 

    # Correlation between noise fields
    if 'noiseCorr' in inputs:
        noiseCorr = np.float(inputs['noiseCorr'])
    else:
        noiseCorr = None

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


    # Is there some signal we want both images to have?
    if 'muBoth' in inputs:

        # Get specification for mu overlap
        muSpecBoth = inputs['muBoth']

        # Reformat center and fwhm if needed
        if 'center' in muSpecBoth:
            muSpecBoth['center']=eval(muSpecBoth['center'])
        if 'fwhm' in muSpecBoth:
            muSpecBoth['fwhm']=eval(muSpecBoth['fwhm'])

    else:

        # Otherwise set it to none
        muSpecBoth = None

    # ID for the configuration
    cfgId = inputs['cfgId']

    # Bootstrap mode
    mode = inputs['mode']

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
    boot_dim = np.array([nSub, 1]) # : COULD MAKE  np.array([batchBoot, nSub, 1])


    # Initialise array for recording the whether set violations occured, for a derived
    # from bootstrapping the true boundary and estimated boundary, respectively. The
    # result in these arrays are based on voxelwise assessment of set condition 
    # violations.
    trueBdry_success1 = np.zeros((nReals,nPvals))
    estBdry_success1 = np.zeros((nReals,nPvals))
    trueBdry_success2 = np.zeros((nReals,nPvals))
    estBdry_success2 = np.zeros((nReals,nPvals))
    trueBdry_success = np.zeros((nReals,nPvals))
    estBdry_success = np.zeros((nReals,nPvals))

    # Initialise array for recording the whether set violations occured, for a derived
    # from bootstrapping the true boundary and estimated boundary, respectively. The
    # result in these arrays are based on interpolation based assessment of set
    # condition violations.
    trueBdry_success_intrp1 = np.zeros((nReals,nPvals))
    estBdry_success_intrp1 = np.zeros((nReals,nPvals))
    trueBdry_success_intrp2 = np.zeros((nReals,nPvals))
    estBdry_success_intrp2 = np.zeros((nReals,nPvals))
    trueBdry_success_intrp = np.zeros((nReals,nPvals))
    estBdry_success_intrp = np.zeros((nReals,nPvals))


    # Initialize for saving times
    times = np.zeros((nReals,1))

    # Loop through realizations
    for r in np.arange(nReals):

        print('r: ', r)
        # -------------------------------------------------------------------
        # Data generation
        # -------------------------------------------------------------------

        # Obtain data
        data1, data2, mu1, mu2 = get_data(muSpec1,muSpec2,noiseSpec1,noiseSpec2,data_dim,noiseCorr,muSpecBoth)

        #print('data shapes: ', data1.shape, data2.shape, mu1.shape, mu2.shape)

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
        # Residuals along dAc1Hat and dAc2Hat
        # -------------------------------------------------------------------

        # Obtain residuals
        resid1 = (data1-muHat1)/sigma1
        resid2 = (data2-muHat2)/sigma2

        # Residuals along Ac boundary
        resid1_dAc1_concat = get_bdry_values_concat(resid1, Ac1_bdry_locs)
        resid2_dAc2_concat = get_bdry_values_concat(resid2, Ac2_bdry_locs)

        # Residuals along AcHat boundary
        resid1_dAcHat1_concat = get_bdry_values_concat(resid1, AcHat1_bdry_locs)
        resid2_dAcHat2_concat = get_bdry_values_concat(resid2, AcHat2_bdry_locs)

        # Delete residuals as they are no longer needed
        #del data1, data2

        # -------------------------------------------------------------------
        # Mu along Ac1 and Ac2 
        # -------------------------------------------------------------------

        # Obtain Mu along Ac
        mu1_dAc1_concat = get_bdry_values_concat(mu1, Ac1_bdry_locs)
        mu2_dAc2_concat = get_bdry_values_concat(mu2, Ac2_bdry_locs)

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


        # print('shapes here : ', mu1_Fc_bdry_concat.shape,resid1_dFc_concat.shape)

        # # Get locations where outer mu1 and mu2 are greater than c
        # dAc1_loc = np.where(mu1_dAc1_concat[0,:,1]>c)[0]
        # dAc2_loc = np.where(mu2_dAc2_concat[0,:,1]>c)[0]

        # print(resid1_dFc_concat[:,d1Fc_loc,:].shape)
        # print(resid2_dFc_concat[:,d2Fc_loc,:].shape)

        # Obtain MuHat along Ac
        muHat1_dAcHat1_concat = get_bdry_values_concat(muHat1, AcHat1_bdry_locs)
        muHat2_dAcHat2_concat = get_bdry_values_concat(muHat2, AcHat2_bdry_locs)

        # print('shapes here : ', muHat1_FcHat_bdry_concat.shape,resid1_FcHat_bdry_concat.shape)

        # # Get locations where outer muHat1 and muHat2 are greater than c
        # dAcHat1_loc = np.where(muHat1_dAcHat1_concat[0,:,1]>c)[0]
        # dAcHat2_loc = np.where(muHat2_dAcHat2_concat[0,:,1]>c)[0]

        # print(resid1_FcHat_bdry_concat[:,d1FcHat_loc,:].shape)
        # # print(resid2_FcHat_bdry_concat[:,d2FcHat_loc,:].shape)

        # # -------------------------------------------------------------------
        # # Residuals along dkFc and dkFcHat boundaries (Array version)
        # # -------------------------------------------------------------------
        # # In this simulation we are bootstrapping the residuals for field 1
        # # along dAc1 intersect Ac2 i.e. along dF where mu2 > c. And vice versa
        # # for field 2.
        # resid1_dAc1_concat = resid1_dAc1_concat[:,dAc1_loc,:]
        # resid2_dAc2_concat = resid2_dAc2_concat[:,dAc2_loc,:]

        # # In this simulation we are bootstrapping the residuals for field 1
        # # along dAcHat1 intersect AcHat2 i.e. along dF where muHat2 > c. And 
        # # vice versa for field 2.
        # resid1_dAcHat1_concat = resid1_dAcHat1_concat[:,dAcHat1_loc,:]
        # resid2_dAcHat2_concat = resid2_dAcHat2_concat[:,dAcHat2_loc,:]

        # # -------------------------------------------------------------------
        # # Mu and MuHat along dAc and dAcHat boundaries (Array version)
        # # -------------------------------------------------------------------
        # mu1_dAc1_concat = mu1_dAc1_concat[:,dAc1_loc,:]
        # mu2_dAc2_concat = mu2_dAc2_concat[:,dAc2_loc,:]

        # # In this simulation we are bootstrapping the residuals for field 1
        # # along dAcHat1 intersect AcHat2 i.e. along dF where muHat2 > c. And 
        # # vice versa for field 2.
        # muHat1_dAcHat1_concat = muHat1_dAcHat1_concat[:,dAcHat1_loc,:]
        # muHat2_dAcHat2_concat = muHat2_dAcHat2_concat[:,dAcHat2_loc,:]

        # -------------------------------------------------------------------
        # Interpolation weights dAc and dAcHat boundaries (Array version)
        # -------------------------------------------------------------------
        dAc1_weights_concat = get_bdry_weights_concat(mu1_dAc1_concat, c)
        dAc2_weights_concat = get_bdry_weights_concat(mu2_dAc2_concat, c)

        dAcHat1_weights_concat = get_bdry_weights_concat(muHat1_dAcHat1_concat, c)
        dAcHat2_weights_concat = get_bdry_weights_concat(muHat2_dAcHat2_concat, c)

        # -------------------------------------------------------------------
        # True and estimated excursion sets
        # -------------------------------------------------------------------
        # Obtain Ac1 and Ac2
        Ac1 = mu1 > c
        Ac2 = mu2 > c

        # Obtain AcHat2 and AcHat2 
        AcHat1 = muHat1 > c
        AcHat2 = muHat2 > c

        # -------------------------------------------------------------------
        # Bootstrap 
        # -------------------------------------------------------------------
        # Initialize empty bootstrap stores
        max_g1_dAc1 = np.zeros(nBoot)
        max_g2_dAc2 = np.zeros(nBoot)
        max_g1_dAcHat1 = np.zeros(nBoot)
        max_g2_dAcHat2 = np.zeros(nBoot)

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

            # Bootstrap residuals along dAc1 and dAc2
            boot_resid1_dAc1_concat = boot_vars*resid1_dAc1_concat
            boot_resid2_dAc2_concat = boot_vars*resid2_dAc2_concat

            # Bootstrap residuals along dAcHat1 and dAcHat2
            boot_resid1_dAcHat1_concat = boot_vars*resid1_dAcHat1_concat
            boot_resid2_dAcHat2_concat = boot_vars*resid2_dAcHat2_concat

            # ------------------------------------------------------------------------------------------
            # True boundary
            # ------------------------------------------------------------------------------------------

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of dAc1. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g1_dAc1_concat = np.zeros(boot_resid1_dAc1_concat.shape[-2:])
            boot_g1_dAc1_concat[...,0] = np.sum(boot_resid1_dAc1_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g1_dAc1_concat[...,1] = np.sum(boot_resid1_dAc1_concat[...,1], axis=0)/np.sqrt(nSub)


            # Sum across subjects to get the bootstrapped a values along
            # the boundary of dAc2. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g2_dAc2_concat = np.zeros(boot_resid2_dAc2_concat.shape[-2:])
            boot_g2_dAc2_concat[...,0] = np.sum(boot_resid2_dAc2_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g2_dAc2_concat[...,1] = np.sum(boot_resid2_dAc2_concat[...,1], axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along dAc1. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma1_boot_dAc1_concat = np.zeros(boot_resid1_dAc1_concat.shape[-2:])
            sigma1_boot_dAc1_concat[...,0] = np.std(boot_resid1_dAc1_concat[...,0], axis=0, ddof=1)
            sigma1_boot_dAc1_concat[...,1] = np.std(boot_resid1_dAc1_concat[...,1], axis=0, ddof=1)

            # Obtain bootstrap standard deviations along dAc2. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma2_boot_dAc2_concat = np.zeros(boot_resid2_dAc2_concat.shape[-2:])
            sigma2_boot_dAc2_concat[...,0] = np.std(boot_resid2_dAc2_concat[...,0], axis=0, ddof=1)
            sigma2_boot_dAc2_concat[...,1] = np.std(boot_resid2_dAc2_concat[...,1], axis=0, ddof=1)

            # Divide by the boostrap standard deviation on dAc1
            boot_g1_dAc1_concat = boot_g1_dAc1_concat/sigma1_boot_dAc1_concat

            # Divide by the boostrap standard deviation on dAc2
            boot_g2_dAc2_concat = boot_g2_dAc2_concat/sigma2_boot_dAc2_concat

            # ------------------------------------------------------------------------------------------
            # Estimated boundary
            # ------------------------------------------------------------------------------------------

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of dAcHat1. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g1_dAcHat1_concat = np.zeros(boot_resid1_dAcHat1_concat.shape[-2:])
            boot_g1_dAcHat1_concat[...,0] = np.sum(boot_resid1_dAcHat1_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g1_dAcHat1_concat[...,1] = np.sum(boot_resid1_dAcHat1_concat[...,1], axis=0)/np.sqrt(nSub)

            # Sum across subjects to get the bootstrapped a values along
            # the boundary of dAcHat2. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. 
            # I am still looking into why this is)
            boot_g2_dAcHat2_concat = np.zeros(boot_resid2_dAcHat2_concat.shape[-2:])
            boot_g2_dAcHat2_concat[...,0] = np.sum(boot_resid2_dAcHat2_concat[...,0], axis=0)/np.sqrt(nSub)
            boot_g2_dAcHat2_concat[...,1] = np.sum(boot_resid2_dAcHat2_concat[...,1], axis=0)/np.sqrt(nSub)

            # Obtain bootstrap standard deviations along dAcHat1. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma1_boot_dAcHat1_concat = np.zeros(boot_resid1_dAcHat1_concat.shape[-2:])
            sigma1_boot_dAcHat1_concat[...,0] = np.std(boot_resid1_dAcHat1_concat[...,0], axis=0, ddof=1)
            sigma1_boot_dAcHat1_concat[...,1] = np.std(boot_resid1_dAcHat1_concat[...,1], axis=0, ddof=1)


            # Obtain bootstrap standard deviations along dAcHat2. (Note: For some reason this is 
            # much faster if performed seperately for each of the last rows. I am still looking
            # into why this is)
            sigma2_boot_dAcHat2_concat = np.zeros(boot_resid2_dAcHat2_concat.shape[-2:])
            sigma2_boot_dAcHat2_concat[...,0] = np.std(boot_resid2_dAcHat2_concat[...,0], axis=0, ddof=1)
            sigma2_boot_dAcHat2_concat[...,1] = np.std(boot_resid2_dAcHat2_concat[...,1], axis=0, ddof=1)

            # Divide by the boostrap standard deviation on dAcHat1
            boot_g1_dAcHat1_concat = boot_g1_dAcHat1_concat/sigma1_boot_dAcHat1_concat

            # Divide by the boostrap standard deviation on dAcHat2
            boot_g2_dAcHat2_concat = boot_g2_dAcHat2_concat/sigma2_boot_dAcHat2_concat

            # ------------------------------------------------------------------------------------------
            # Interpolation
            # ------------------------------------------------------------------------------------------

            # Interpolation for dAc1 and dAc2 boundary
            boot_g1_dAc1_concat = get_bdry_vals_interpolated_concat(boot_g1_dAc1_concat,dAc1_weights_concat)
            boot_g2_dAc2_concat = get_bdry_vals_interpolated_concat(boot_g2_dAc2_concat,dAc2_weights_concat)

            # Interpolation for dAcHat1 and dAcHat2 boundary
            boot_g1_dAcHat1_concat = get_bdry_vals_interpolated_concat(boot_g1_dAcHat1_concat,dAcHat1_weights_concat)
            boot_g2_dAcHat2_concat = get_bdry_vals_interpolated_concat(boot_g2_dAcHat2_concat,dAcHat2_weights_concat)


            # Get the maximum along dAc1 if it exists
            if np.prod(boot_g1_dAc1_concat.shape) > 0:
                max_g1_dAc1[b] = np.max(np.abs(boot_g1_dAc1_concat)) 

            # Get the maximum along dAc2 if it exists
            if np.prod(boot_g2_dAc2_concat.shape) > 0:
                max_g2_dAc2[b] = np.max(np.abs(boot_g2_dAc2_concat)) 

            # Get the maximum along dAcHat1 if it exists
            if np.prod(boot_g1_dAcHat1_concat.shape) > 0:
                max_g1_dAcHat1[b] = np.max(np.abs(boot_g1_dAcHat1_concat)) 


            # Get the maximum along dAcHat2 if it exists
            if np.prod(boot_g2_dAcHat2_concat.shape) > 0:
                max_g2_dAcHat2[b] = np.max(np.abs(boot_g2_dAcHat2_concat)) 


            # -------------------------------------------------------------------------------------

            
        t2 = time.time()
        # print('Bootstrap time: ', t2-t1)

        # -------------------------------------------------------------------
        # Obtaining a from percentiles of the max distribution
        # -------------------------------------------------------------------

        # Drop the instances where the boundary length was zero
        max_g1_dAc1 = max_g1_dAc1[max_g1_dAc1!=0]
        max_g2_dAc2 = max_g2_dAc2[max_g2_dAc2!=0]
        max_g1_dAcHat1 = max_g1_dAcHat1[max_g1_dAcHat1!=0]
        max_g2_dAcHat2 = max_g2_dAcHat2[max_g2_dAcHat2!=0]

        # If we have recorded values get their quantiles
        if (np.prod(max_g1_dAc1.shape) > 0):
            # Get the a estimates for the true boundary
            a_trueBdry1 = np.percentile(max_g1_dAc1, 100*p).reshape(nPvals,1,1,1)
        else:
            # Set to inf by default
            a_trueBdry1 = np.Inf*np.ones((nPvals,1,1,1))

        # If we have recorded values get their quantiles
        if (np.prod(max_g2_dAc2.shape) > 0):
            # Get the a estimates for the true boundary
            a_trueBdry2 = np.percentile(max_g2_dAc2, 100*p).reshape(nPvals,1,1,1)
        else:
            # Set to inf by default
            a_trueBdry2 = np.Inf*np.ones((nPvals,1,1,1))

        # If we have recorded values get their quantiles
        if (np.prod(max_g1_dAcHat1.shape) > 0):
            # Get the a estimates for the estimated boundary
            a_estBdry1 = np.percentile(max_g1_dAcHat1, 100*p).reshape(nPvals,1,1,1) # [pvals, 1, [1 for _ in dim>1]]
        else:
            # Set to inf by default
            a_estBdry1 = np.Inf*np.ones((nPvals,1,1,1))

        # If we have recorded values get their quantiles
        if (np.prod(max_g2_dAcHat2.shape) > 0):
            # Get the a estimates for the estimated boundary
            a_estBdry2 = np.percentile(max_g2_dAcHat2, 100*p).reshape(nPvals,1,1,1) # [pvals, 1, [1 for _ in dim>1]]
        else:
            # Set to inf by default
            a_estBdry2 = np.Inf*np.ones((nPvals,1,1,1))

        # Reformat them to an array form useful for boolean operation
        a_trueBdry1 = np.concatenate((-a_trueBdry1,a_trueBdry1),axis=1)
        a_trueBdry2 = np.concatenate((-a_trueBdry2,a_trueBdry2),axis=1)
        a_estBdry1 = np.concatenate((-a_estBdry1,a_estBdry1),axis=1)
        a_estBdry2 = np.concatenate((-a_estBdry2,a_estBdry2),axis=1)

        # print('a')
        # print(a_trueBdry)
        # print(a_estBdry)

        # Get the statistic field which defined Achat^{+/-,1/2}
        g1 = ((muHat1-c)/(sigma1*tau))
        g1 = g1.reshape(g1.shape[-2],g1.shape[-1])
        g2 = ((muHat2-c)/(sigma2*tau))
        g2 = g2.reshape(g2.shape[-2],g2.shape[-1])

        # Obtain AcHat^+ and AcHat^- based on a from the true boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        AcHat1_pm_trueBdry = g1 >= a_trueBdry1
        AcHat2_pm_trueBdry = g2 >= a_trueBdry2

        # Obtain AcHat^+ and AcHat^- based on a from the estimated boundary. This variable
        # has axes corresponding to [pvalue, plus/minus, field dimensions]
        AcHat1_pm_estBdry = g1 >= a_estBdry1
        AcHat2_pm_estBdry = g2 >= a_estBdry2


        # End timer (we have the confidence sets now)
        t2 = time.time()

        # Save time 
        times[r,:] = t2-t1
        
        # -------------------------------------------------------------------
        # Images
        # -------------------------------------------------------------------
        if r==1 and inputs['figGen']:

            # Indices for 0.8,0.9 and 0.95
            pInds = [16,18,19]

            for pInd in pInds:

                # P value of interest
                pVal = p[pInd]

                # AcHat1 plus (based on estimated bdry)
                plt.figure(16)
                plt.imshow(1*AcHat1_pm_estBdry[pInd,0,...])
                plt.savefig(os.path.join(figDir, 'AcHat1_plus_estBdry_p'+str(int(100*pVal))+'_cfg'+str(cfgId)+'.png'))

                # AcHat1 minus (based on estimated bdry)
                plt.figure(17)
                plt.imshow(1*AcHat1_pm_estBdry[pInd,1,...])
                plt.savefig(os.path.join(figDir, 'AcHat1_minus_estBdry_p'+str(int(100*pVal))+'_cfg'+str(cfgId)+'.png'))

                # AcHat1 plus (based on true bdry)
                plt.figure(18)
                plt.imshow(1*AcHat1_pm_trueBdry[pInd,0,...])
                plt.savefig(os.path.join(figDir, 'AcHat1_plus_trueBdry_p'+str(int(100*pVal))+'_cfg'+str(cfgId)+'.png'))

                # AcHat1 minus (based on true bdry)
                plt.figure(19)
                plt.imshow(1*AcHat1_pm_trueBdry[pInd,1,...])
                plt.savefig(os.path.join(figDir, 'AcHat1_minus_trueBdry_p'+str(int(100*pVal))+'_cfg'+str(cfgId)+'.png'))


                # AcHat2 plus (based on estimated bdry)
                plt.figure(20)
                plt.imshow(1*AcHat2_pm_estBdry[pInd,0,...])
                plt.savefig(os.path.join(figDir, 'AcHat2_plus_estBdry_p'+str(int(100*pVal))+'_cfg'+str(cfgId)+'.png'))

                # AcHat2 minus (based on estimated bdry)
                plt.figure(21)
                plt.imshow(1*AcHat2_pm_estBdry[pInd,1,...])
                plt.savefig(os.path.join(figDir, 'AcHat2_minus_estBdry_p'+str(int(100*pVal))+'_cfg'+str(cfgId)+'.png'))

                # AcHat2 plus (based on true bdry)
                plt.figure(22)
                plt.imshow(1*AcHat2_pm_trueBdry[pInd,0,...])
                plt.savefig(os.path.join(figDir, 'AcHat2_plus_trueBdry_p'+str(int(100*pVal))+'_cfg'+str(cfgId)+'.png'))

                # AcHat2 minus (based on true bdry)
                plt.figure(23)
                plt.imshow(1*AcHat2_pm_trueBdry[pInd,1,...])
                plt.savefig(os.path.join(figDir, 'AcHat2_minus_trueBdry_p'+str(int(100*pVal))+'_cfg'+str(cfgId)+'.png'))

        # -------------------------------------------------------------------
        # Some set logic to work out violations
        # -------------------------------------------------------------------

        # Obtain AcHat1^+\AcHat1 based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        AcHat1p_sub_Ac1_trueBdry = AcHat1_pm_trueBdry[:,1,...] & ~Ac1[...]

        # Obtain AcHat1^+\AcHat1 based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        AcHat1p_sub_Ac1_estBdry = AcHat1_pm_estBdry[:,1,...] & ~Ac1[...]

        # Obtain Ac1\AcHat1^- based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        Ac1_sub_AcHat1m_trueBdry = Ac1[...] & ~AcHat1_pm_trueBdry[:,0,...]

        # Obtain Ac1\AcHat1^- based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        Ac1_sub_AcHat1m_estBdry = Ac1[...] & ~AcHat1_pm_estBdry[:,0,...]

        # Obtain AcHat2^+\AcHat2 based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        AcHat2p_sub_Ac2_trueBdry = AcHat2_pm_trueBdry[:,1,...] & ~Ac2[...]

        # Obtain AcHat2^+\AcHat2 based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        AcHat2p_sub_Ac2_estBdry = AcHat2_pm_estBdry[:,1,...] & ~Ac2[...]

        # Obtain Ac2\AcHat2^- based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        Ac2_sub_AcHat2m_trueBdry = Ac2[...] & ~AcHat2_pm_trueBdry[:,0,...]

        # Obtain Ac2\AcHat2^- based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        Ac2_sub_AcHat2m_estBdry = Ac2[...] & ~AcHat2_pm_estBdry[:,0,...]


        # Record if we saw a violation in the true boundary based sets
        trueBdry_success1[r,:] = 1-(np.any(AcHat1p_sub_Ac1_trueBdry,axis=(1,2)) | np.any(Ac1_sub_AcHat1m_trueBdry,axis=(1,2))) # : AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success1[r,:] = 1-(np.any(AcHat1p_sub_Ac1_estBdry,axis=(1,2)) | np.any(Ac1_sub_AcHat1m_estBdry,axis=(1,2)))# : AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the true boundary based sets
        trueBdry_success2[r,:] = 1-(np.any(AcHat2p_sub_Ac2_trueBdry,axis=(1,2)) | np.any(Ac2_sub_AcHat2m_trueBdry,axis=(1,2))) # : AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success2[r,:] = 1-(np.any(AcHat2p_sub_Ac2_estBdry,axis=(1,2)) | np.any(Ac2_sub_AcHat2m_estBdry,axis=(1,2)))# : AXES WONT WORK FOR 3D ATM

        # -------------------------------------------------------------------
        # Intersection images
        # -------------------------------------------------------------------

        # Work out intersection masks
        AcHatCap_pm_trueBdry = AcHat1_pm_trueBdry & AcHat2_pm_trueBdry
        AcHatCap_pm_estBdry = AcHat1_pm_estBdry & AcHat2_pm_estBdry
        AcCap = Ac1 & Ac2

        # -------------------------------------------------------------------
        # Some set logic to work out violations for intersection
        # -------------------------------------------------------------------

        # Obtain AcHat1^+\AcHat1 based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        AcHatCapp_sub_AcCap_trueBdry = AcHatCap_pm_trueBdry[:,1,...] & ~AcCap[...]

        # Obtain AcHat1^+\AcHat1 based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        AcHatCapp_sub_AcCap_estBdry = AcHatCap_pm_estBdry[:,1,...] & ~AcCap[...]

        # Obtain Ac1\AcHat1^- based on the true boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        AcCap_sub_AcHatCapm_trueBdry = AcCap[...] & ~AcHatCap_pm_trueBdry[:,0,...]

        # Obtain Ac1\AcHat1^- based on the estimated boundary. This variable
        # has axes corresponding to [pvalue, field dimensions]
        AcCap_sub_AcHatCapm_estBdry = AcCap[...] & ~AcHatCap_pm_estBdry[:,0,...]

        # Record if we saw a violation in the true boundary based sets
        trueBdry_success[r,:] = 1-(np.any(AcHatCapp_sub_AcCap_trueBdry,axis=(1,2)) | np.any(AcCap_sub_AcHatCapm_trueBdry,axis=(1,2))) # : AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success[r,:] = 1-(np.any(AcHatCapp_sub_AcCap_estBdry,axis=(1,2)) | np.any(AcCap_sub_AcHatCapm_estBdry,axis=(1,2)))# : AXES WONT WORK FOR 3D ATM


        # -------------------------------------------------------------------
        # Get stat along the Ac boundaries
        # -------------------------------------------------------------------


        # Obtain g1 and g2 along the boundary for Ac1 and Ac2
        g1_dAc1_concat = get_bdry_values_concat(g1.reshape(mu1.shape), Ac1_bdry_locs)
        g2_dAc2_concat = get_bdry_values_concat(g2.reshape(mu2.shape), Ac2_bdry_locs)

        # print('g1 shape ', g1.shape)
        # print('g2 shape ', g2.shape)
        # print('Ac1_bdry_locs ', Ac1_bdry_locs)
        # print('Ac2_bdry_locs ', Ac2_bdry_locs)


        # Interpolation
        g1_dAc1_concat = get_bdry_vals_interpolated_concat(g1_dAc1_concat,dAc1_weights_concat)
        g2_dAc2_concat = get_bdry_vals_interpolated_concat(g2_dAc2_concat,dAc2_weights_concat)

        # Reshape a bit
        g1_dAc1_concat = g1_dAc1_concat.reshape(g1_dAc1_concat.shape[-2],g1_dAc1_concat.shape[-1])
        g2_dAc2_concat = g2_dAc2_concat.reshape(g2_dAc2_concat.shape[-2],g2_dAc2_concat.shape[-1])

        # -------------------------------------------------------------------
        # Intersection locations
        # -------------------------------------------------------------------

        # Get boolean maps for the boundary of Fc
        Fc_bdry_maps = get_bdry_maps(np.minimum(mu1,mu2), c)

        # Get coordinates for the boundary of Fc
        Fc_bdry_locs = get_bdry_locs(Fc_bdry_maps)

        # Obtain Mu along Fc (We need this to estimate interpolated regions)
        mu1_Fc_bdry_concat = get_bdry_values_concat(mu1, Fc_bdry_locs)
        mu2_Fc_bdry_concat = get_bdry_values_concat(mu2, Fc_bdry_locs)

        # Obtain g1 and g2 along the boundary for Fc
        g1_dFc_concat = get_bdry_values_concat(g1, Fc_bdry_locs)
        g2_dFc_concat = get_bdry_values_concat(g2, Fc_bdry_locs)
        
        # Get locations where outer mu1 and mu2 are greater than c
        d1Fc_loc = np.where((mu2_Fc_bdry_concat[0,:,1]>c))[0]
        d2Fc_loc = np.where((mu1_Fc_bdry_concat[0,:,1]>c))[0]
        d12Fc_loc = np.where((mu2_Fc_bdry_concat[0,:,1]<= c)*(mu1_Fc_bdry_concat[0,:,1]<=c))[0]

        # print('g1 shape ', g1.shape)
        # print('g2 shape ', g2.shape)
        # print('Ac1_bdry_locs ', Ac1_bdry_locs)
        # print('Ac2_bdry_locs ', Ac2_bdry_locs)

        # -------------------------------------------------------------------
        # Interpolation weights
        # -------------------------------------------------------------------

        # In this simulation we are bootstrapping the residuals for field 1
        # along dAc1 intersect Ac2 i.e. along dF where mu2 > c. And vice versa
        # for field 2.
        mu1_d1Fc_concat = mu1_Fc_bdry_concat[:,d1Fc_loc,:]
        mu2_d2Fc_concat = mu2_Fc_bdry_concat[:,d2Fc_loc,:]

        # For intersect boundary
        mu1_d12Fc_concat = mu1_Fc_bdry_concat[:,d12Fc_loc,:]
        mu2_d12Fc_concat = mu2_Fc_bdry_concat[:,d12Fc_loc,:]

        # Obtain the weights along the boundary for dkFc
        d1Fc_bdry_weights_concat = get_bdry_weights_concat(mu1_d1Fc_concat, c)
        d2Fc_bdry_weights_concat = get_bdry_weights_concat(mu2_d2Fc_concat, c)
        d12Fc_mu1_bdry_weights_concat = get_bdry_weights_concat(mu1_d12Fc_concat, c)
        d12Fc_mu2_bdry_weights_concat = get_bdry_weights_concat(mu2_d12Fc_concat, c)

        # -------------------------------------------------------------------
        # g1 and g2 on boundary
        # -------------------------------------------------------------------

        # g1 and g2 along d1Fc
        g1_d1Fc_concat = g1_dFc_concat[:,d1Fc_loc,:]
        g2_d1Fc_concat = g2_dFc_concat[:,d1Fc_loc,:]

        # g1 and g2 along d2Fc
        g1_d2Fc_concat = g1_dFc_concat[:,d2Fc_loc,:]
        g2_d2Fc_concat = g2_dFc_concat[:,d2Fc_loc,:]

        # g1 and g2 along d12Fc
        g1_d12Fc_concat = g1_dFc_concat[:,d12Fc_loc,:]
        g2_d12Fc_concat = g2_dFc_concat[:,d12Fc_loc,:]

        # Interpolation for d1Fc boundary (we perform this interpolation based on the weights of
        # mu1 as we are on the boundary of Ac1)
        g1_d1Fc_concat = get_bdry_vals_interpolated_concat(g1_d1Fc_concat,d1Fc_bdry_weights_concat)
        g2_d1Fc_concat = get_bdry_vals_interpolated_concat(g2_d1Fc_concat,d1Fc_bdry_weights_concat)

        # Interpolation for d2Fc boundary (we perform this interpolation based on the weights of
        # mu2 as we are on the boundary of Ac2)
        g1_d2Fc_concat = get_bdry_vals_interpolated_concat(g1_d2Fc_concat,d2Fc_bdry_weights_concat)
        g2_d2Fc_concat = get_bdry_vals_interpolated_concat(g2_d2Fc_concat,d2Fc_bdry_weights_concat)


        # Interpolation for d1Fc and d2Fc boundary
        g1_d12Fc_concat = get_bdry_vals_interpolated_concat(g1_d12Fc_concat,d12Fc_mu1_bdry_weights_concat)
        g2_d12Fc_concat = get_bdry_vals_interpolated_concat(g2_d12Fc_concat,d12Fc_mu2_bdry_weights_concat)

        # Reshape a bit
        g1_d1Fc_concat = g1_d1Fc_concat.reshape(g1_d1Fc_concat.shape[-2], g1_d1Fc_concat.shape[-1])
        g2_d1Fc_concat = g2_d1Fc_concat.reshape(g2_d1Fc_concat.shape[-2], g2_d1Fc_concat.shape[-1])
        g1_d2Fc_concat = g1_d2Fc_concat.reshape(g1_d2Fc_concat.shape[-2], g1_d2Fc_concat.shape[-1])
        g2_d2Fc_concat = g2_d2Fc_concat.reshape(g2_d2Fc_concat.shape[-2], g2_d2Fc_concat.shape[-1])
        g1_d12Fc_concat = g1_d12Fc_concat.reshape(g1_d12Fc_concat.shape[-2], g1_d12Fc_concat.shape[-1])
        g2_d12Fc_concat = g2_d12Fc_concat.reshape(g2_d12Fc_concat.shape[-2], g2_d12Fc_concat.shape[-1])


        # Get g1 and g2 along dFc
        g1_FcBdry = np.concatenate((g1_d1Fc_concat,g1_d2Fc_concat,g1_d12Fc_concat),axis=-1)
        g1_FcBdry = g1_FcBdry.reshape(g1_FcBdry.shape[-2],g1_FcBdry.shape[-1])

        g2_FcBdry = np.concatenate((g2_d1Fc_concat,g2_d2Fc_concat,g2_d12Fc_concat),axis=-1)
        g2_FcBdry = g2_FcBdry.reshape(g2_FcBdry.shape[-2],g2_FcBdry.shape[-1])

        # -------------------------------------------------------------------
        # Check whether there were any boundary violations using interpolated
        # boundary values (checking if voxels had values corresponding to no
        # violations, etc) (For intersection)
        # -------------------------------------------------------------------

        # Perform lower check on stat map using thresholds based on the
        # estimated boundary
        bdry_lowerCheck_estBdry_cap1 = g1_dFc_concat >= a_estBdry1[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # estimated boundary
        bdry_upperCheck_estBdry_cap1 = g1_dFc_concat <= a_estBdry1[:,1,:,0]

        # Perform lower check on stat map using thresholds based on the
        # true boundary
        bdry_lowerCheck_trueBdry_cap1 = g1_dFc_concat >= a_trueBdry1[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # true boundary
        bdry_upperCheck_trueBdry_cap1 = g1_dFc_concat <= a_trueBdry1[:,1,:,0]

        # Perform lower check on stat map using thresholds based on the
        # estimated boundary
        bdry_lowerCheck_estBdry_cap2 = g2_dFc_concat >= a_estBdry2[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # estimated boundary
        bdry_upperCheck_estBdry_cap2 = g2_dFc_concat <= a_estBdry2[:,1,:,0]

        # Perform lower check on stat map using thresholds based on the
        # true boundary
        bdry_lowerCheck_trueBdry_cap2 = g2_dFc_concat >= a_trueBdry2[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # true boundary
        bdry_upperCheck_trueBdry_cap2 = g2_dFc_concat <= a_trueBdry2[:,1,:,0]


        # Success in both instances
        trueBdry_success_intrp[r,:] = np.all(bdry_lowerCheck_trueBdry_cap1,axis=(1)) & np.all(bdry_upperCheck_trueBdry_cap1,axis=(1)) & \
                                      np.all(bdry_lowerCheck_trueBdry_cap2,axis=(1)) & np.all(bdry_upperCheck_trueBdry_cap2,axis=(1))
        estBdry_success_intrp[r,:] = np.all(bdry_lowerCheck_estBdry_cap1,axis=(1)) & np.all(bdry_upperCheck_estBdry_cap1,axis=(1)) & \
                                     np.all(bdry_lowerCheck_estBdry_cap2,axis=(1)) & np.all(bdry_upperCheck_estBdry_cap2,axis=(1))

        # -------------------------------------------------------------------
        # Check whether there were any boundary violations using interpolated
        # boundary values (checking if voxels had values corresponding to no
        # violations, etc) (For Ac1 and Ac2 seperately)
        # -------------------------------------------------------------------

        # Perform lower check on stat map using thresholds based on the
        # estimated boundary
        bdry_lowerCheck_estBdry1 = g1_dAc1_concat >= a_estBdry1[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # estimated boundary
        bdry_upperCheck_estBdry1 = g1_dAc1_concat <= a_estBdry1[:,1,:,0]

        # Perform lower check on stat map using thresholds based on the
        # true boundary
        bdry_lowerCheck_trueBdry1 = g1_dAc1_concat >= a_trueBdry1[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # true boundary
        bdry_upperCheck_trueBdry1 = g1_dAc1_concat <= a_trueBdry1[:,1,:,0]

        # Perform lower check on stat map using thresholds based on the
        # estimated boundary
        bdry_lowerCheck_estBdry2 = g2_dAc2_concat >= a_estBdry2[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # estimated boundary
        bdry_upperCheck_estBdry2 = g2_dAc2_concat <= a_estBdry2[:,1,:,0]

        # Perform lower check on stat map using thresholds based on the
        # true boundary
        bdry_lowerCheck_trueBdry2 = g2_dAc2_concat >= a_trueBdry2[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # true boundary
        bdry_upperCheck_trueBdry2 = g2_dAc2_concat <= a_trueBdry2[:,1,:,0]

        # -------------------------------------------------------------------
        # Work out whether simulation observed successful sets.
        # -------------------------------------------------------------------
        # Record if we saw a violation in the true boundary based sets
        trueBdry_success_intrp1[r,:] = np.all(bdry_lowerCheck_trueBdry1,axis=(1)) & np.all(bdry_upperCheck_trueBdry1,axis=(1))# : AXES WONT WORK FOR 3D ATM
        trueBdry_success_intrp2[r,:] = np.all(bdry_lowerCheck_trueBdry2,axis=(1)) & np.all(bdry_upperCheck_trueBdry2,axis=(1))# : AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success_intrp1[r,:] = np.all(bdry_lowerCheck_estBdry1,axis=(1)) & np.all(bdry_upperCheck_estBdry1,axis=(1)) # : AXES WONT WORK FOR 3D ATM
        estBdry_success_intrp2[r,:] = np.all(bdry_lowerCheck_estBdry2,axis=(1)) & np.all(bdry_upperCheck_estBdry2,axis=(1)) # : AXES WONT WORK FOR 3D ATM


    # For the interpolated boundary success checks, we still need to do the 
    # voxelwise checks as well. This will take care of that.
    trueBdry_success_intrp1 = trueBdry_success_intrp1*trueBdry_success1
    estBdry_success_intrp1 = estBdry_success_intrp1*estBdry_success1

    # For the interpolated boundary success checks, we still need to do the 
    # voxelwise checks as well. This will take care of that.
    trueBdry_success_intrp2 = trueBdry_success_intrp2*trueBdry_success2
    estBdry_success_intrp2 = estBdry_success_intrp2*estBdry_success2

    # For the interpolated boundary success checks, we still need to do the 
    # voxelwise checks as well. This will take care of that.
    trueBdry_success_intrp = trueBdry_success_intrp*trueBdry_success
    estBdry_success_intrp = estBdry_success_intrp*estBdry_success

    # Coverage probabilities
    coverage_trueBdry1 = np.mean(trueBdry_success1,axis=0)
    coverage_estBdry1 = np.mean(estBdry_success1,axis=0)

    # Coverage probabilities
    coverage_trueBdry2 = np.mean(trueBdry_success2,axis=0)
    coverage_estBdry2 = np.mean(estBdry_success2,axis=0)

    # Coverage probabilities
    coverage_trueBdry = np.mean(trueBdry_success,axis=0)
    coverage_estBdry = np.mean(estBdry_success,axis=0)

    # Coverage probabilities
    coverage_trueBdry_intrp1 = np.mean(trueBdry_success_intrp1,axis=0)
    coverage_estBdry_intrp1 = np.mean(estBdry_success_intrp1,axis=0)

    # Coverage probabilities
    coverage_trueBdry_intrp2 = np.mean(trueBdry_success_intrp2,axis=0)
    coverage_estBdry_intrp2 = np.mean(estBdry_success_intrp2,axis=0)

    # Coverage probabilities
    coverage_trueBdry_intrp = np.mean(trueBdry_success_intrp,axis=0)
    coverage_estBdry_intrp = np.mean(estBdry_success_intrp,axis=0)

    # Make results folder
    if not os.path.exists(os.path.join(simDir, 'RawResults')):
        os.mkdir(os.path.join(simDir, 'RawResults'))

    # Save the violations to a file
    append_to_file(os.path.join(simDir, 'RawResults', 'trueSuccess1.csv'), trueBdry_success1) # Successes based on the true boundary (assessed without interpolation)
    append_to_file(os.path.join(simDir, 'RawResults', 'estSuccess1.csv'), estBdry_success1) # Successes based on the interpolated boundary (assessed without interpolation) 
    append_to_file(os.path.join(simDir, 'RawResults', 'trueSuccess_intrp1.csv'), trueBdry_success_intrp1) # Successes based on the true boundary (assessed with interpolation)
    append_to_file(os.path.join(simDir, 'RawResults', 'estSuccess_intrp1.csv'), estBdry_success_intrp1) # Successes based on the interpolated boundary (assessed with interpolation) 

    append_to_file(os.path.join(simDir, 'RawResults', 'trueSuccess2.csv'), trueBdry_success2) # Successes based on the true boundary (assessed without interpolation)
    append_to_file(os.path.join(simDir, 'RawResults', 'estSuccess2.csv'), estBdry_success2) # Successes based on the interpolated boundary (assessed without interpolation) 
    append_to_file(os.path.join(simDir, 'RawResults', 'trueSuccess_intrp2.csv'), trueBdry_success_intrp2) # Successes based on the true boundary (assessed with interpolation)
    append_to_file(os.path.join(simDir, 'RawResults', 'estSuccess_intrp2.csv'), estBdry_success_intrp2) # Successes based on the interpolated boundary (assessed with interpolation) 

    append_to_file(os.path.join(simDir, 'RawResults', 'trueSuccess.csv'), trueBdry_success) # Successes based on the true boundary (assessed without interpolation)
    append_to_file(os.path.join(simDir, 'RawResults', 'estSuccess.csv'), estBdry_success) # Successes based on the interpolated boundary (assessed without interpolation) 
    append_to_file(os.path.join(simDir, 'RawResults', 'trueSuccess_intrp.csv'), trueBdry_success_intrp) # Successes based on the true boundary (assessed with interpolation)
    append_to_file(os.path.join(simDir, 'RawResults', 'estSuccess_intrp.csv'), estBdry_success_intrp) # Successes based on the interpolated boundary (assessed with interpolation) 

    append_to_file(os.path.join(simDir, 'RawResults', 'times.csv'), times) # Times for bootstrap

    # Save the computation times
    t2overall = time.time()
    append_to_file(os.path.join(simDir, 'RawResults', 'computationTime.csv'), np.array([t2overall-t1overall]))

#SpatialSims_2mu('/home/tommaullin/Documents/ConfRes/tmp/sim15/sim15/cfgs/cfg578.yml')
