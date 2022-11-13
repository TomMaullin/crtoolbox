import os
import time
import numpy as np
from lib.generateData import *
from lib.boundary import *
from lib.fileio import *
import yaml
import matplotlib.pyplot as plt

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

    # Make results folder
    if not os.path.exists(os.path.join(simDir, 'RawResults')):
        os.mkdir(os.path.join(simDir, 'RawResults'))

    # Save the violations to a file
    append_to_file(os.path.join(simDir, 'RawResults', 'trueSuccess.csv'), trueBdry_success) # Successes based on the true boundary (assessed without interpolation)
    append_to_file(os.path.join(simDir, 'RawResults', 'estSuccess.csv'), estBdry_success) # Successes based on the interpolated boundary (assessed without interpolation) 
    append_to_file(os.path.join(simDir, 'RawResults', 'trueSuccess_intrp.csv'), trueBdry_success_intrp) # Successes based on the true boundary (assessed with interpolation)
    append_to_file(os.path.join(simDir, 'RawResults', 'estSuccess_intrp.csv'), estBdry_success_intrp) # Successes based on the interpolated boundary (assessed with interpolation) 
    append_to_file(os.path.join(simDir, 'RawResults', 'times.csv'), times) # Times for bootstrap
