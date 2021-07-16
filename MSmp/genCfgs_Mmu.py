import os
import numpy as np
import yaml

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
    inputs['nReals'] = 500

    # Add number of bootstraps
    inputs['nBoot'] = 5000

    # Number of fields, m
    ms = np.arange(2,6)

    # This determines how close together the circles are.
    # (spacing r is the distance from each circle center to the
    # center of the image)
    spacing_r = 25

    # These are our sample sizes:
    nSubs = np.linspace(40,500,24)
    fg_nSubs = np.array([100,300,500]) 

    # Add mode
    inputs['mode'] = 3

    # ==========================================================================
    #
    # Simulation 1: Circles in a circle
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in m circles equally spread 
    # in a circle, as we vary m
    #
    # ==========================================================================
    if simNo==1:
    
        # Add threshold c
        inputs['c'] = 2

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

    # ==========================================================================
    #
    # Simulation 2: Circles in a circle (lower FWHM)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in m circles equally spread 
    # in a circle, as we vary m
    #
    # ==========================================================================
    if simNo==2:
    
        # Add threshold c
        inputs['c'] = 2/3
        
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
                mus['mu'+str(i+1)]['mag'] = 1

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