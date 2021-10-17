import os
import time
import numpy as np
from lib.generateData import *
from lib.boundary import *
from lib.fileio import *
import yaml
from scipy.ndimage.measurements import label
import matplotlib.pyplot as plt

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

m = 5
alphas=list(powerset(np.arange(m)))


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

def get_data_1field(muSpec,noiseSpec,dim):

    # Obtain the noise fields
    noise = get_noise(noiseSpec, dim)

    # Obtain mu
    mu = get_mu(muSpec, dim)
    
    # Create the data
    data = mu + noise

    # Return the data and mu
    return(data,mu)

def SpatialSims_Mmu(ipath):

    # -----------------------------------------------------------------------
    # Load in inputs 
    # -----------------------------------------------------------------------
    # Read in file
    with open(ipath, 'r') as stream:
        inputs = yaml.load(stream,Loader=yaml.FullLoader)

    # Get output directory
    OutDir = inputs['OutDir']

    # ID for the configuration
    cfgId = inputs['cfgId']

    # Get number of subjects
    nSub = int(inputs['nSub'])

    # Get simulation number
    simNo = int(inputs['simNo'])

    # Work out how many fields we have
    m = np.int(inputs['m'])

    # Get number of subjects
    nSub = int(inputs['nSub'])

    # Get number of simulation realizations
    nReals = int(inputs['nReals'])

    # Get number of bootstraps
    nBoot = int(inputs['nBoot'])

    # Get Threshold
    c = np.float(inputs['c'])

    # Get p values
    p = eval(inputs['p'])

    # Get the number of p-values we're looking at
    nPvals = len(p)

    # Get number of bootstraps
    if 'tau' in inputs:
        tau = eval(inputs['tau'])
    else:
        tau = 1/np.sqrt(nSub)

    # Dimensions of simulated data
    data_dim = np.array([nSub, 100,100])

    # Dimensions of bootstrap variables
    boot_dim = np.array([nSub, 1]) 

    # Initialise array for recording the whether set violations occured, for a derived
    # from bootstrapping the true boundary and estimated boundary, respectively. The
    # result in these arrays are based on voxelwise assessment of set condition 
    # violations.
    trueBdry_success = np.zeros((nReals,nPvals))
    estBdry_success = np.zeros((nReals,nPvals))

    # Simulation directory
    simDir = os.path.join(OutDir, 'sim'+str(simNo), 'cfg' + str(cfgId))
    if not os.path.exists(simDir):
        os.mkdir(simDir)

    # ----------------------------------------------------------------
    # Mu specification
    # ----------------------------------------------------------------

    # Empty dict to hold mu specifications
    muSpec = {}

    # Loop through values of i
    for i in np.arange(m):

        # Get specification for mu i
        muSpeci = inputs['mus']['mu'+str(i+1)]

        # Reformat center and fwhm if needed
        if 'center' in muSpeci:
            muSpeci['center']=eval(muSpeci['center'])
        if 'fwhm' in muSpeci:
            muSpeci['fwhm']=eval(muSpeci['fwhm'])

        # Save muSpec
        muSpec[str(i+1)] = muSpeci

    # ----------------------------------------------------------------
    # Noise specification
    # ----------------------------------------------------------------

    # Empty dict to hold noise specifications
    noiseSpec = {}

    # Loop through values of i
    for i in np.arange(m):

        # Get noise spec
        noiseSpeci = inputs['noises']['noise'+str(i+1)]

        # Reformat FWHM for the noise
        noiseSpeci['FWHM'] = str2vec(noiseSpeci['FWHM']) 

        # Save noiseSpec
        noiseSpec[str(i+1)] = noiseSpeci

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

        print('Realization: ', r)

        # Make a structure to hold the true and estimated boundary locations
        true_bdry_locs = {}
        est_bdry_locs = {}

        # Make a structure to hold the true and estimated boundary weights
        true_bdry_weights = {}
        est_bdry_weights = {}

        # Make a structure to hold the true and estimated boundary weights in
        # array form (it's current handy to have both forms close to hand but
        # this will likely be changed in later versions of the code)
        true_bdry_weights_concat = {}
        est_bdry_weights_concat = {}

        for i in np.arange(m):

            # ----------------------------------------------------------------
            # Data generation
            # ----------------------------------------------------------------

            # Obtain data
            data, mu = get_data_1field(muSpec[str(i+1)],noiseSpec[str(i+1)],data_dim)

            # Save mus
            if i == 0:
                mus = np.array(mu)
            else:
                mus = np.concatenate((mus,mu),axis=0)

            # Combine data
            if i == 0:
                datas = np.array(data.reshape(1,*(data.shape)))
            else:
                datas = np.concatenate((datas,data.reshape(1,*(data.shape))),axis=0)

            # plt.figure(int(i))
            # plt.imshow(1*mu[0,:,:])


            # -------------------------------------------------------------------
            # Mean and standard deviation estimates
            # -------------------------------------------------------------------

            # Obtain mu estimate
            muHat = np.mean(data, axis=0).reshape(mu.shape)

            # Save muHats
            if i == 0:
                muHats = np.array(muHat)
            else:
                muHats = np.concatenate((muHats,muHat),axis=0)

            # Obtain sigma
            sigma = np.std(data, axis=0).reshape(mu.shape)

            # Save sigmas
            if i == 0:
                sigmas = np.array(sigma)
            else:
                sigmas = np.concatenate((sigmas,sigma),axis=0)

            # -------------------------------------------------------------------
            # Boundary locations for Aci
            # -------------------------------------------------------------------
            # Get boolean map for the boundary of Aci
            Ac_bdry_map = get_bdry_maps(mu, c)

            # Get coordinates for the boundary of Aci
            Ac_bdry_locs = get_bdry_locs(Ac_bdry_map)

            # Save boundary locations
            true_bdry_locs['Ac'+str(i+1)] = Ac_bdry_locs

            # Delete map as we no longer need it
            del Ac_bdry_map

            # -------------------------------------------------------------------
            # Boundary locations for AcHati
            # -------------------------------------------------------------------
            # Get boolean maps for the boundary of AcHat
            AcHat_bdry_map = get_bdry_maps(muHat, c)

            # Get coordinates for the boundary of AcHat
            AcHat_bdry_locs = get_bdry_locs(AcHat_bdry_map)

            # Save boundary locations
            est_bdry_locs['AcHat'+str(i+1)] = AcHat_bdry_locs

            # Delete map as we no longer need it
            del AcHat_bdry_map

            # -------------------------------------------------------------------
            # Interpolation weights for Aci boundary (Dict version)
            # -------------------------------------------------------------------
            # Obtain the values along the boundary for Aci
            Ac_bdry_vals = get_bdry_values(mu, Ac_bdry_locs)

            # Obtain the weights along the boundary for Aci
            Ac_bdry_weights = get_bdry_weights(Ac_bdry_vals, c)

            # Save boundary weights
            true_bdry_weights['Ac'+str(i+1)] = Ac_bdry_weights

            # Delete values as we no longer need them
            del Ac_bdry_vals

            # -------------------------------------------------------------------
            # Interpolation weights for Aci boundary (Array version)
            # -------------------------------------------------------------------
            # Obtain the values along the boundary for Aci
            Ac_bdry_vals_concat = get_bdry_values_concat(mu, Ac_bdry_locs)

            # Obtain the weights along the boundary for Aci
            Ac_bdry_weights_concat = get_bdry_weights_concat(Ac_bdry_vals_concat, c)

            # Save boundary weights
            true_bdry_weights_concat['Ac'+str(i+1)] = Ac_bdry_weights_concat

            # Delete values as we no longer need them
            del Ac_bdry_vals_concat

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
        # Get minimum fields
        # -------------------------------------------------------------------
        # This is named cap as the excursion set of the minimum field is
        # the intersection of all fields (\cap in latex)
        cap_mu = np.amin(mus,axis=0)
        cap_muHat = np.amin(muHats,axis=0)

        # -------------------------------------------------------------------
        # Boundary locations for Fc 
        # -------------------------------------------------------------------
        # Get boolean map for the boundary of Fc
        Fc_bdry_map = get_bdry_maps(cap_mu, c)

        # Get coordinates for the boundary of Fc
        Fc_bdry_locs = get_bdry_locs(Fc_bdry_map)

        # Save boundary locations
        true_bdry_locs['Fc'] = Fc_bdry_locs

        # Delete maps as we no longer need them
        del Fc_bdry_map

        # -------------------------------------------------------------------
        # Boundary locations for FcHat
        # -------------------------------------------------------------------
        # Get boolean map for the boundary of FcHat
        FcHat_bdry_map = get_bdry_maps(cap_muHat, c)

        # Get coordinates for the boundary of FcHat
        FcHat_bdry_locs = get_bdry_locs(FcHat_bdry_map)

        # Save boundary locations
        est_bdry_locs['FcHat'] = FcHat_bdry_locs

        # Delete maps as we no longer need them
        del FcHat_bdry_map

        # -------------------------------------------------------------------
        # Interpolation weights for Fc boundary (Dict version)
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for Fc
        Fc_bdry_vals = get_bdry_values(cap_mu, Fc_bdry_locs)

        # Obtain the weights along the boundary for Fc
        Fc_bdry_weights = get_bdry_weights(Fc_bdry_vals, c)

        # Save boundary weights
        true_bdry_weights['Fc'] = Fc_bdry_weights

        # Delete values as we no longer need them
        del Fc_bdry_vals

        # -------------------------------------------------------------------
        # Interpolation weights for Fc boundary (Array version)
        # -------------------------------------------------------------------
        # Obtain the values along the boundary for Fc
        Fc_bdry_vals_concat = get_bdry_values_concat(cap_mu, Fc_bdry_locs)

        # Obtain the weights along the boundary for Fc
        Fc_bdry_weights_concat = get_bdry_weights_concat(Fc_bdry_vals_concat, c)

        # Save boundary weights
        true_bdry_weights_concat['Fc'] = Fc_bdry_weights_concat

        # Delete values as we no longer need them
        del Fc_bdry_vals_concat

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

        # Empty dicts to store residuals
        resids_dFc = {}
        resids_dFcHat = {}

        # Empty dicts to store mu and muHat
        mu_dFc = {}
        muHat_dFcHat = {}

        # Loop through to get residuals
        for i in np.arange(m):

            # -------------------------------------------------------------------
            # Residuals along dFc and dFcHat
            # -------------------------------------------------------------------

            # Obtain residuals
            resid = (datas[i,...]-muHats[i,...])/sigmas[i,...] # NEED DATA TO BE SAVED FIRST

            # Residuals along Fc boundary
            resid_dFc_concat = get_bdry_values_concat(resid, Fc_bdry_locs)

            # Save residuals
            resids_dFc['field'+str(i+1)] = resid_dFc_concat

            # Residuals along FcHat boundary
            resid_dFcHat_concat = get_bdry_values_concat(resid, FcHat_bdry_locs)

            # Save residuals
            resids_dFcHat['field'+str(i+1)] = resid_dFcHat_concat

            # -------------------------------------------------------------------
            # Mu along dFc
            # -------------------------------------------------------------------

            # Obtain Mu along Fc
            mu_Fc_bdry_concat = get_bdry_values_concat(mus[i,...], Fc_bdry_locs)

            # Save mu
            mu_dFc['field'+str(i+1)] = mu_Fc_bdry_concat

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
        dalphaFc_locs = {}
        dalphaFcHat_locs = {}

        # Loop through alphas
        for alpha in alphas:

            # Loop through possible elements of alpha
            for i in (np.arange(m)+1):

                # Get indices representing where mui and muhati are less
                # than or equal to c if i is in alpha.
                if i in alpha:

                    # Mu
                    in_dalphaFc_i = (mu_dFc['field'+str(i)][:,1] <= c)

                    # Muhat
                    in_dalphaFcHat_i = (muHat_dFcHat['field'+str(i)][:,1] <= c)

                # Get indices representing where mui is greater than c
                # if i is not in alpha.
                else:

                    # Mu
                    in_dalphaFc_i = (mu_dFc['field'+str(i)][:,1] > c)

                    # Muhat
                    in_dalphaFcHat_i = (muHat_dFcHat['field'+str(i)][:,1] > c)

                # print('alpha: ', alpha)
                # print('i: ', i)
                # print(in_dalpha_i)

                # Get a boolean index telling us whether each location in dFc
                # belongs to d^\alpha Fc.
                if i == 1:

                    # Initial running product for mu
                    in_dalphaFc = np.array(in_dalphaFc_i)

                    # Initial running product for muHat
                    in_dalphaFcHat = np.array(in_dalphaFcHat_i)

                else:

                    # Update running product for mu
                    in_dalphaFc = in_dalphaFc*in_dalphaFc_i

                    # Update running product for muhat
                    in_dalphaFcHat = in_dalphaFcHat*in_dalphaFcHat_i

            # Convert boolean 1/0s into coordinates
            dalphaFc_loc = np.where(in_dalphaFc)[0]
            dalphaFcHat_loc = np.where(in_dalphaFcHat)[0]

            # Save locations
            dalphaFc_locs[np.array2string(alpha)] = dalphaFc_loc
            dalphaFcHat_locs[np.array2string(alpha)] = dalphaFcHat_loc

        # -------------------------------------------------------------------
        # Get residuals, mu and muhat along boundary partitions
        # -------------------------------------------------------------------

        # Empty dicts for residuals, mu and muHat
        resids_dFc_partitioned = {}
        resids_dFcHat_partitioned = {}
        mu_dFc_partitioned = {}
        muHat_dFcHat_partitioned = {}

        # Loop through boundary partitions
        for alpha in alphas:

            # New empty dicts for Fc
            resids_dalphaFc = {}
            mu_dalphaFc = {}

            # New empty dicts for FcHat
            resids_dalphaFcHat = {}
            muHat_dalphaFcHat = {}

            # Get dalpha locations
            dalphaFc_loc = dalphaFc_locs[np.array2string(alpha)]
            dalphaFcHat_loc = dalphaFcHat_locs[np.array2string(alpha)]

            # Loop through i in alpha getting values for interpolation
            for i in alpha:

                # ------------------------------------------------------
                # Residuals and mu on true boundary; dFc
                # ------------------------------------------------------

                # Get residuals for field i along dFc
                residsi_dFc = resids_dFc['field'+str(i)]

                # Save residuals for field i along dalphaFc
                resids_dalphaFc[str(i)] = residsi_dFc[:,dalphaFc_loc,:]

                # Get mu for field i along dFc
                mui_dFc = mu_dFc['field'+str(i)]

                # Save residuals for field i along dalphaFc
                mu_dalphaFc[str(i)] = mui_dFc[dalphaFc_loc,:]

                # ------------------------------------------------------
                # Residuals and mu on estimates boundary; dFcHat
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
            resids_dFc_partitioned[np.array2string(alpha)] = resids_dalphaFc
            resids_dFcHat_partitioned[np.array2string(alpha)] = resids_dalphaFcHat

            # Save mu and muhat for alpha
            mu_dFc_partitioned[np.array2string(alpha)] = mu_dalphaFc
            muHat_dFcHat_partitioned[np.array2string(alpha)] = muHat_dalphaFcHat

        # -------------------------------------------------------------------
        # Get weights from mu and muhat along boundary partitions for 
        # interpolation
        # -------------------------------------------------------------------

        # Empty dicts for weights
        weights_dFc = {}
        weights_dFcHat = {}

        # Loop through boundary partitions
        for alpha in alphas:

            # New empty dicts for weights
            weights_dalphaFc = {}
            weights_dalphaFcHat = {}

            # Loop through i in alpha getting values for interpolation
            for i in alpha:

                # Get mui and muhati along dalpha Fc
                mui_dalphaFc = mu_dFc_partitioned[np.array2string(alpha)][str(i)]
                muHati_dalphaFcHat = muHat_dFcHat_partitioned[np.array2string(alpha)][str(i)]

                # For true boundary
                dalphaFc_mui_weights = get_bdry_weights_concat(mui_dalphaFc, c)

                # For estimated boundary
                dalphaFcHat_muHati_weights = get_bdry_weights_concat(muHati_dalphaFcHat, c)

                # Save weights
                weights_dalphaFc[str(i)] = dalphaFc_mui_weights
                weights_dalphaFcHat[str(i)] = dalphaFcHat_muHati_weights

            # Save weights
            weights_dFc[np.array2string(alpha)] = weights_dalphaFc
            weights_dFcHat[np.array2string(alpha)] = weights_dalphaFcHat


        # -------------------------------------------------------------------
        # True and estimated excursion sets
        # -------------------------------------------------------------------
        # Obtain Ac for all i
        Ac = mus > c

        # Obtain Fc
        Fc = cap_mu > c

        # Obtain AcHat for all i
        AcHat = muHats > c

        # Obtain FcHat
        FcHat = cap_muHat > c


        # -------------------------------------------------------------------
        # Bootstrap 
        # -------------------------------------------------------------------
        # Initialize empty bootstrap stores for supremum along each boundary
        # e.g. supg_dFc[str(alpha)][str(i)] will give the supremum of g^i 
        #      along the boundary d\alpha F_c
        supg_dFc = {}
        supg_dFcHat = {}

        # Loop through boundary partitions
        for alpha in alphas:

            # Initialize
            supg_dalphaFc = {}
            supg_dalphaFcHat = {}

            # Loop through i in alpha getting values for interpolation
            for i in alpha:

                # Save empty bootstrap stores
                supg_dalphaFc[str(i)] = np.zeros(nBoot)
                supg_dalphaFcHat[str(i)] = np.zeros(nBoot)

            # Save empty bootstrap stores
            supg_dFc[np.array2string(alpha)] = supg_dalphaFc
            supg_dFcHat[np.array2string(alpha)] = supg_dalphaFcHat

        # Initalize empty bootstrap stores for the min
        min_supg_dFc = {}
        min_supg_dFcHat = {}

        # Loop through boundary partitions
        for alpha in alphas:

            # Save empty bootstrap stores
            min_supg_dFc[np.array2string(alpha)] = np.zeros(nBoot)
            min_supg_dFcHat[np.array2string(alpha)] = np.zeros(nBoot)

        # Save empty bootstrap stores
        min_supg_dFc['max'] = np.zeros(nBoot)
        min_supg_dFcHat['max'] = np.zeros(nBoot)

        t1 = time.time()
        # For each bootstrap record the max of the residuals along the
        # boundary
        for b in np.arange(nBoot):

            # print('Bootstrap: ', b)

            # Obtain bootstrap variables
            boot_vars = 2*np.random.randint(0,2,boot_dim)-1

            # Reshape for broadcasting purposes (extra axis refers to the fact we have
            # inner and outer boundary values in the last axes of resid_FsdGc_bdry_concat
            # and resid_FsdGcHat_bdry_concat)
            boot_vars = boot_vars.reshape((*boot_vars.shape),1)
            
            # -------------------------------------------------------------------------
            # Bootstrap residuals along dalphaFc
            # -------------------------------------------------------------------------
            boot_resids_dFc = {}
            boot_resids_dFcHat = {}

            # Loop through boundary partitions
            for alpha in alphas:

                # New empty dicts for Fc
                boot_resids_dalphaFc = {}

                # New empty dicts for FcHat
                boot_resids_dalphaFcHat = {}

                # New empty dict for gi along dalphaFc
                boot_g_dalphaFc = {}

                # New empty dict for gi along dalphaFcHat
                boot_g_dalphaFcHat = {}

                # Loop through i in alpha getting gi interpolated
                for i in alpha:

                    # ------------------------------------------------------
                    # Get residuals true and estimated boundary
                    # ------------------------------------------------------

                    # Get original residuals of mui along boundary Aci
                    residsi_dalphaFc = resids_dFc_partitioned[np.array2string(alpha)][str(i)]
                    residsi_dalphaFcHat = resids_dFcHat_partitioned[np.array2string(alpha)][str(i)]

                    # ------------------------------------------------------
                    # Bootstrap residuals
                    # ------------------------------------------------------

                    # Multiply by rademacher variables
                    boot_residsi_dalphaFc = boot_vars*residsi_dalphaFc
                    boot_residsi_dalphaFcHat = boot_vars*residsi_dalphaFcHat

                    # ------------------------------------------------------
                    # Get gi along dalpha Fc
                    # ------------------------------------------------------

                    # Sum across subjects to get the bootstrapped a values along
                    # the boundary of dalphaFc. (Note: For some reason this is 
                    # much faster if performed seperately for each of the last rows. 
                    # I am still looking into why this is)
                    boot_gi_dalphaFc = np.zeros(boot_residsi_dalphaFc.shape[-2:])
                    boot_gi_dalphaFc[...,0] = np.sum(boot_residsi_dalphaFc[...,0], axis=0)/np.sqrt(nSub)
                    boot_gi_dalphaFc[...,1] = np.sum(boot_residsi_dalphaFc[...,1], axis=0)/np.sqrt(nSub)

                    # Obtain bootstrap standard deviations for mu i along dalphaFc. 
                    # (Note: For some reason this is much faster if performed seperately
                    # for each of the last rows. I am still looking into why this is)
                    boot_sigmai_dalphaFc = np.zeros(boot_residsi_dalphaFc.shape[-2:])
                    boot_sigmai_dalphaFc[...,0] = np.std(boot_residsi_dalphaFc[...,0], axis=0, ddof=1)
                    boot_sigmai_dalphaFc[...,1] = np.std(boot_residsi_dalphaFc[...,1], axis=0, ddof=1)

                    # Divide by the boostrap standard deviation of mui
                    boot_gi_dalphaFc = boot_gi_dalphaFc/boot_sigmai_dalphaFc

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
                    # Interpolate along dalpha Fc
                    # ------------------------------------------------------

                    # Get weights
                    dalphaFc_mui_bdry_weights = weights_dFc[np.array2string(alpha)][str(i)]

                    # Interpolation for gi along dalphaFc
                    boot_gi_dalphaFc = get_bdry_vals_interpolated_concat(boot_gi_dalphaFc,dalphaFc_mui_bdry_weights)

                    # ------------------------------------------------------
                    # Interpolate along dalpha FcHat
                    # ------------------------------------------------------

                    # Get weights
                    dalphaFcHat_muHati_bdry_weights = weights_dFcHat[np.array2string(alpha)][str(i)]

                    # Interpolation for gi along dalphaFc
                    boot_gi_dalphaFcHat = get_bdry_vals_interpolated_concat(boot_gi_dalphaFcHat,dalphaFcHat_muHati_bdry_weights)

                    # ------------------------------------------------------
                    # Get supremum along dalphaFc
                    # ------------------------------------------------------

                    # Get the maximum of gi along dalphaFc if it exists
                    if np.prod(boot_gi_dalphaFc.shape) > 0:

                        boot_g_dalphaFc[str(i)] = boot_gi_dalphaFc

                    # ------------------------------------------------------
                    # Get supremum along dalphaFcHat
                    # ------------------------------------------------------

                    # Get the maximum of gi along dalphaFcHat if it exists
                    if np.prod(boot_gi_dalphaFcHat.shape) > 0:

                        boot_g_dalphaFcHat[str(i)] = boot_gi_dalphaFcHat

                # Boolean telling us if this is the first field we are looking at
                first = True

                # -----------------------------------------------------------------
                # Get sup_{dalpha Fc} |min_{alpha} g^i|
                # -----------------------------------------------------------------

                # Loop through i in alpha getting min(gi)
                for i in alpha:

                    # Check if we have boundary values for gi
                    if str(i) in boot_g_dalphaFc:

                        # If this is the first time initalize the min(g) array
                        if first:

                            # Initialize elementwise minimum across gi of boot_g
                            boot_ming_dalphaFc = boot_g_dalphaFc[str(i)] 

                            # We are no longer looking at the first
                            first = False

                        else:

                            # Get elementwise minimum across gi of boot_g
                            boot_ming_dalphaFc = np.minimum(boot_ming_dalphaFc,boot_g_dalphaFc[str(i)])

                # If we saw some boundary values (i.e. we saw at least the first i),
                # then the array was initialized
                if not first:

                    # Save sup(|min(gi)|) along true dalpha Fc
                    boot_sup_ming_dalphaFc = np.max(np.abs(boot_ming_dalphaFc))

                    # Save bootstrap result
                    min_supg_dFc[np.array2string(alpha)][b] = boot_sup_ming_dalphaFc


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

                # Check if we have boundary values for dalpha Fc
                if np.array2string(alpha) in min_supg_dFc:

                    # Update the minimum value we've seen
                    min_supg_dFc['max'][b] = np.maximum(min_supg_dFc['max'][b],min_supg_dFc[np.array2string(alpha)][b])

                # Check if we have boundary values for dalpha FcHat
                if np.array2string(alpha) in min_supg_dFcHat:

                    # Update the minimum value we've seen
                    min_supg_dFcHat['max'][b] = np.maximum(min_supg_dFcHat['max'][b],min_supg_dFcHat[np.array2string(alpha)][b])

        # End timer (we have the confidence sets now)
        t2 = time.time()

        # -------------------------------------------------------------------
        # Obtaining a from percentiles of the max distribution
        # -------------------------------------------------------------------

        # Drop the instances where the boundary length was zero
        min_supg_dFc['max'] = min_supg_dFc['max'][min_supg_dFc['max']!=0]
        min_supg_dFcHat['max'] = min_supg_dFcHat['max'][min_supg_dFcHat['max']!=0]


        # If we have recorded values get their quantiles
        if (np.prod(min_supg_dFc['max'].shape) > 0):
            # Get the a estimates for the true boundary
            a_trueBdry = np.percentile(min_supg_dFc['max'], 100*p).reshape(nPvals,1,1,1)
        else:
            # Set to inf by default
            a_trueBdry = np.Inf*np.ones((nPvals,1,1,1))

        # If we have recorded values get their quantiles
        if (np.prod(min_supg_dFcHat['max'].shape) > 0):
            # Get the a estimates for the estimated boundary
            a_estBdry = np.percentile(min_supg_dFcHat['max'], 100*p).reshape(nPvals,1,1,1) # [pvals, 1, [1 for _ in dim>1]]
        else:
            # Set to inf by default
            a_estBdry = np.Inf*np.ones((nPvals,1,1,1))

        # Reformat them to an array form useful for boolean operation
        a_trueBdry = np.concatenate((-a_trueBdry,a_trueBdry),axis=1)
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

        # Save time 
        times[r,:] = t2-t1

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

        #print('anys: ', np.any(FcHatp_sub_Fc_trueBdry,axis=(1,2)), np.any(Fc_sub_FcHatm_trueBdry,axis=(1,2)))
        # Record if we saw a violation in the true boundary based sets
        trueBdry_success[r,:] = 1-(np.any(FcHatp_sub_Fc_trueBdry,axis=(1,2)) | np.any(Fc_sub_FcHatm_trueBdry,axis=(1,2))) # : AXES WONT WORK FOR 3D ATM

        # Record if we saw a violation in the estimated boundary based sets
        estBdry_success[r,:] = 1-(np.any(FcHatp_sub_Fc_estBdry,axis=(1,2)) | np.any(Fc_sub_FcHatm_estBdry,axis=(1,2)))# : AXES WONT WORK FOR 3D ATM

        # -------------------------------------------------------------------
        # Get stat along the Fc boundary
        # -------------------------------------------------------------------

        # Empty dict to store g along true boundary
        g_dFc = {}

        # Loop through fields
        for i in np.arange(m):

            # Get the values for gi along dFc
            g_dFc[str(i+1)] = get_bdry_values_concat(g[i,...], Fc_bdry_locs)

        # Empty dict to store interpolate g along true boundary
        g_dFc_interp = {}

        # Boolean to tell us if this is the first alpha we've looked at
        Firstalpha = True

        # Interpolate g
        for alpha in alphas:

            # Initial structure to hold interpolated g
            g_dalphaFc_interp = {}

            # Boolean to tell us if this is the first value of i we're looking at
            Firsti = True

            # Check if we have recorded values for this boundary
            if np.array2string(alpha) in dalphaFc_locs:

                # Loop through all fields relevent to boundary.
                for i in alpha:

                    # Get gi on dFc
                    gi_dFc = g_dFc[str(i)] 

                    # Get locations for dalphaFc
                    dalphaFc_loc = dalphaFc_locs[np.array2string(alpha)]

                    # Get gi on dalphaFc
                    gi_dalphaFc = gi_dFc[dalphaFc_loc,:]

                    # Get weights for interpolation
                    dalphaFc_weightsi = weights_dFc[np.array2string(alpha)][str(i)]

                    # Interpolate gi
                    gi_dalphaFc_interp = get_bdry_vals_interpolated_concat(gi_dalphaFc,dalphaFc_weightsi)

                    # If this is the first i we've looked at in alpha
                    if Firsti:

                        # set the minimum interpolated values to gi_dalphaFc_interp
                        ming_dalphaFc_interp = gi_dalphaFc_interp 

                        # Set first to false
                        Firsti = False

                    else:

                        # Take the minimum between current gi_dalphaFc_interp and previous
                        ming_dalphaFc_interp = np.minimum(ming_dalphaFc_interp, gi_dalphaFc_interp) 

                # -------------------------------------------------------------------
                # Get stat along the Fc boundary
                # -------------------------------------------------------------------

                # If this is the first alpha we've looked at, set the stat array to ming_dalphaFc_interp
                if Firstalpha:

                    stat_dFc = ming_dalphaFc_interp

                    # Now we've seen an alpha
                    Firstalpha = False

                # Otherwise concatenate it on
                else:

                    stat_dFc = np.concatenate((ming_dalphaFc_interp,stat_dFc),axis=-1)

        # -------------------------------------------------------------------
        # Check whether there were any boundary violations using interpolated
        # boundary values (checking if voxels had values corresponding to no
        # violations, etc)
        # -------------------------------------------------------------------

        # Perform lower check on stat map using thresholds based on the
        # estimated boundary
        bdry_lowerCheck_estBdry = stat_dFc >= a_estBdry[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # estimated boundary
        bdry_upperCheck_estBdry = stat_dFc <= a_estBdry[:,1,:,0]

        # Perform lower check on stat map using thresholds based on the
        # true boundary
        bdry_lowerCheck_trueBdry = stat_dFc >= a_trueBdry[:,0,:,0]

        # Perform upper check on stat map using thresholds based on the
        # true boundary
        bdry_upperCheck_trueBdry = stat_dFc <= a_trueBdry[:,1,:,0]


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
    
    print('coverage_estBdry_intrp: ', coverage_estBdry_intrp)

    # Make results folder
    if not os.path.exists(os.path.join(simDir, 'RawResults')):
        os.mkdir(os.path.join(simDir, 'RawResults'))

    # Save the violations to a file
    append_to_file(os.path.join(simDir, 'RawResults', 'trueSuccess.csv'), trueBdry_success) # Successes based on the true boundary (assessed without interpolation)
    append_to_file(os.path.join(simDir, 'RawResults', 'estSuccess.csv'), estBdry_success) # Successes based on the interpolated boundary (assessed without interpolation) 
    append_to_file(os.path.join(simDir, 'RawResults', 'trueSuccess_intrp.csv'), trueBdry_success_intrp) # Successes based on the true boundary (assessed with interpolation)
    append_to_file(os.path.join(simDir, 'RawResults', 'estSuccess_intrp.csv'), estBdry_success_intrp) # Successes based on the interpolated boundary (assessed with interpolation) 
    append_to_file(os.path.join(simDir, 'RawResults', 'times.csv'), times) # Times for bootstrap
