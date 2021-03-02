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


        print('shapes here : ', mu1_Fc_bdry_concat.shape,resid1_dFc_concat.shape)

        # Get locations where outer mu1 and mu2 are greater than c
        d1Fc_loc = np.where(mu2_Fc_bdry_concat[0,:,1]>c)[0]
        d2Fc_loc = np.where(mu1_Fc_bdry_concat[0,:,1]>c)[0]
        d12Fc_loc = np.where((mu2_Fc_bdry_concat[0,:,1]<= c)*(mu1_Fc_bdry_concat[0,:,1]<=c))[0]

        print(resid1_dFc_concat[:,d1Fc_loc,:].shape)
        print(resid2_dFc_concat[:,d2Fc_loc,:].shape)

        # Obtain MuHat along Fc
        muHat1_FcHat_bdry_concat = get_bdry_values_concat(muHat1, FcHat_bdry_locs)
        muHat2_FcHat_bdry_concat = get_bdry_values_concat(muHat2, FcHat_bdry_locs)

        print('shapes here : ', muHat1_FcHat_bdry_concat.shape,resid1_FcHat_bdry_concat.shape)

        # Get locations where outer muHat1 and muHat2 are greater than c
        d1FcHat_loc = np.where(muHat2_FcHat_bdry_concat[0,:,1]>c)[0]
        d2FcHat_loc = np.where(muHat1_FcHat_bdry_concat[0,:,1]>c)[0]
        d12FcHat_loc = np.where((muHat2_FcHat_bdry_concat[0,:,1]<= c)*(muHat1_FcHat_bdry_concat[0,:,1]<=c))[0]

        print(resid1_FcHat_bdry_concat[:,d1FcHat_loc,:].shape)
        print(resid2_FcHat_bdry_concat[:,d2FcHat_loc,:].shape)

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


        print('marker')
        print(mu1_d1Fc_concat.shape)
        print(mu1_Fc_bdry_concat.shape)
        print(d2Fc_loc.shape)
        print(mu2_d2Fc_concat.shape)
        print(mu2_Fc_bdry_concat.shape)
        print(d1Fc_loc.shape)
        print(muHat1_d1FcHat_concat.shape)

        print('marker2')
        print(d1FcHat_loc.shape)
        print(d2FcHat_loc.shape)
        print(d12FcHat_loc.shape)
        print(muHat1_FcHat_bdry_concat.shape)

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

                # Get the maximum along d12Fc if it exists
                if np.prod(boot_g1_d12Fc_concat.shape) > 0:

                    # Minimum of both g1 and g2 on d12Fc boundary (note: the absolute values must be outside
                    # the minimum in this mode, unlike in mode 2 where they were inside)
                    min_g1g2_d12Fc[b] = np.max(np.abs(np.minimum(boot_g1_d12Fc_concat,boot_g2_d12Fc_concat)))

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
            max_g_dFc = np.maximum(max_g1_d1Fc,max_g2_d2Fc,max_g1g2_d12Fc)

            # Get the maximum needed for the estimated boundary
            max_g_dFcHat = np.maximum(max_g1_d1FcHat,max_g2_d2FcHat,max_g1g2_d12FcHat)


        # In mode 3, we bootstrap G1 along d1Fc and G2 along d2Fc and min(G1,G2) along d12Fc
        if mode == 3:

            # Get the maximum needed for the true boundary
            max_g_dFc = np.maximum(max_g1_d1Fc,max_g2_d2Fc,min_g1g2_d12Fc)

            # Get the maximum needed for the estimated boundary
            max_g_dFcHat = np.maximum(max_g1_d1FcHat,max_g2_d2FcHat,min_g1g2_d12FcHat)

        t2 = time.time()
        print('Bootstrap time: ', t2-t1)

        # -------------------------------------------------------------------
        # Obtaining a from percentiles of the max distribution
        # -------------------------------------------------------------------

        # Get the a estimates for the true boundaries and estimated boundaries
        a_trueBdry = np.percentile(max_g_dFc, 100*p).reshape(nPvals,1,1,1)
        a_estBdry = np.percentile(max_g_dFcHat, 100*p).reshape(nPvals,1,1,1) # [pvals, 1, [1 for _ in dim>1]]

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

        #print(Fc_sub_FcHatm_estBdry.shape, Fc.shape, FcHat_pm_estBdry.shape)

        # -------------------------------------------------------------------
        # Check whether there were any boundary violations using voxelwise
        # set logic (checking if voxels existed in one set but not another, 
        # etc)
        # -------------------------------------------------------------------
        # Record if we saw a violation in the true boundary based sets
        trueBdry_success[r,:] = 1-(np.any(FcHatp_sub_Fc_trueBdry,axis=(1,2)) | np.any(Fc_sub_FcHatm_trueBdry,axis=(1,2))) # : AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success[r,:] = 1-(np.any(FcHatp_sub_Fc_estBdry,axis=(1,2)) | np.any(Fc_sub_FcHatm_estBdry,axis=(1,2)))# : AXES WONT WORK FOR 3D ATM

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

    print('Coverage: ', coverage_estBdry_intrp)
    print('Coverage: ', coverage_trueBdry_intrp)

    # Make results folder
    if not os.path.exists(os.path.join(simDir, 'RawResults')):
        os.mkdir(os.path.join(simDir, 'RawResults'))

    # Save the violations to a file
    append_to_file(os.path.join(simDir, 'RawResults', 'trueSuccess.csv'), trueBdry_success) # Successes based on the true boundary (assessed without interpolation)
    append_to_file(os.path.join(simDir, 'RawResults', 'estSuccess.csv'), estBdry_success) # Successes based on the interpolated boundary (assessed without interpolation) 
    append_to_file(os.path.join(simDir, 'RawResults', 'trueSuccess_intrp.csv'), trueBdry_success_intrp) # Successes based on the true boundary (assessed with interpolation)
    append_to_file(os.path.join(simDir, 'RawResults', 'estSuccess_intrp.csv'), estBdry_success_intrp) # Successes based on the interpolated boundary (assessed with interpolation) 

    # Save the computation times
    t2overall = time.time()
    append_to_file(os.path.join(simDir, 'RawResults', 'computationTime.csv'), np.array([t2overall-t1overall]))

#SpatialSims_2mu('/home/tommaullin/Documents/ConfSets/sim1/cfgs/cfg623.yml')