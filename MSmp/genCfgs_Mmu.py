import os
import numpy as np
import yaml

def generateCfgs(OutDir, simNo):

    # Make simulation directory
    simDir = os.path.join(OutDir, 'sim'+str(simNo))
    if not os.path.exists(simDir):
        os.mkdir(simDir)

    # Make directory to store configuration files
    if not os.path.exists(os.path.join(simDir,'cfgs')):
        os.mkdir(os.path.join(simDir,'cfgs'))

    # ======================================================================================
    # Input generation
    # ======================================================================================

    # Make simulation directory
    simDir = os.path.join(OutDir, 'sim'+str(simNo))
    if not os.path.exists(simDir):
        os.mkdir(simDir)

    # Make directory to store configuration files
    if not os.path.exists(os.path.join(simDir,'cfgs')):
        os.mkdir(os.path.join(simDir,'cfgs'))

    # -------------------------------------------------------------------------------------
    # Basic inputs
    # -------------------------------------------------------------------------------------

    # New empty inputs structure
    inputs={}

    # Add output directory
    inputs['OutDir'] = OutDir

    # Add simulation number (helps identify simulation)
    inputs['simNo'] = simNo

    # Add range of p-values
    inputs['p'] = 'np.linspace(0,1,21)'

    # Add tau expression
    inputs['tau'] = '1/np.sqrt(nSub)'

    # Add number of realizations
    inputs['nReals'] = 50#0

    # Add number of bootstraps
    inputs['nBoot'] = 500#0

    # Number of fields, m
    ms = np.arange(2,8)

    # This determines how close together the circles are.
    # (spacing r is the distance from each circle center to the
    # center of the image)
    spacing_r = 25

    # These are our sample sizes:
    nSubs = np.linspace(40,500,24)
    fg_nSubs = np.array([100,300,500]) 

    # Add threshold c
    inputs['c'] = 2

    # Add mode
    inputs['mode'] = 3

    # Id for config file
    cfgId = 1

    # Loop through all m values
    for m in ms:

        # Add m
        inputs['m']=int(m)

        # ---------------------------------------------------------------
        # Mus
        # ---------------------------------------------------------------
        # Create empty specifications
        mus = {}

        # Loop through mus, adding each field in turn
        for i in np.arange(m):

            # New empty dict
            mus['mu'+str(i+1)]={}

            # Mu type
            mus['mu'+str(i+1)]['type'] = 'circle2D' 

            # Mu FWHM
            mus['mu'+str(i+1)]['fwhm'] = 'np.array([5,5])'

            # Mu r
            mus['mu'+str(i+1)]['r'] = 40

            # Mu magnitude
            mus['mu'+str(i+1)]['mag'] = 3

            # Get some evenly spaced center points
            centers = circle_points(np.array([spacing_r]),np.array([m]))

            # Mu center
            mus['mu'+str(i+1)]['center'] = 'np.'+repr(centers[i,:].astype(np.int))

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
            noises['noise'+str(i+1)]['FWHM'] = '[0, 3, 3]'

            # Add type
            noises['noise'+str(i+1)]['type'] = 'homogen'

        # Save mus and noises
        inputs['mus']=mus
        inputs['noises']=noises

        # Loop through all nSub settings
        for nSub in nSubs:

            # Add nSub to inputs
            inputs['nSub'] = int(nSub)

            # Save cfg ID (handy to have around)
            inputs['cfgId'] = int(cfgId)

            # Record if we want to save figures for this design or not
            if (nSub in fg_nSubs):

                # In this case we do want to save  figures
                inputs['figGen']=1

            else:

                # In this case we do want to save  figures
                inputs['figGen']=0

            # Save the yml
            with open(os.path.join(simDir,'cfgs','cfg'+str(cfgId)+'.yml'), 'w') as outfile:
                yaml.dump(inputs, outfile, default_flow_style=False)

            # Incremement cfgID
            cfgId = cfgId + 1


    # Delete fields which vary across simulation
    del inputs['mus'], inputs['noises'], inputs['cfgId'], inputs['nSub']

    #--------------------------------------------------------------------------------------
    # Save the baseline configuration (Note: This must be the last file output as the bash
    # script used for running simulations on the cluster will take this files existence as
    # a sign to run the next stage of the simulations)
    #--------------------------------------------------------------------------------------
    # Save the yml
    with open(os.path.join(simDir,'cfgs','baseline_cfg.yml'), 'w') as outfile:
        yaml.dump(inputs, outfile, default_flow_style=False)