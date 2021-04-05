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


    # ==========================================================================
    # Basic inputs for all simulations
    # ==========================================================================

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

    # Create mu1 and mu2 specification
    mu1 = {}
    mu2 = {}

    # Create noise1 and noise2 specification
    noise1 = {}
    noise2 = {}

    # These are our sample sizes:
    nSubs = np.linspace(40,500,24)
    fg_nSubs = np.array([100,300,500])  

    # ==========================================================================
    #
    # Simulation 1: Circles moving closer
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two circles of 
    # equal diameter close to one another. For this reason, we vary the circles
    # center and, as usual, the number of subjects.
    #
    # ==========================================================================
    if simNo==1:

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # Add threshold c
        inputs['c'] = 2

        # Add mode
        inputs['mode'] = 1

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'circle2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center (we only vary mu2)
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'circle2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for center in centers:

            # Add mu2 center
            mu2['center']= 'np.array(['+str(center)+',0])'

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

    # ==========================================================================
    #
    # Simulation 2: Circles moving closer (but lower fwhm)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two circles of 
    # equal diameter close to one another. For this reason, we vary the circles
    # center and, as usual, the number of subjects. This differs from the 
    # previous setting in terms of smoothing and fwhm.
    #
    # ==========================================================================
    if simNo==2:

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # Add threshold c
        inputs['c'] = 2/3

        # Add mode
        inputs['mode'] = 1

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'circle2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 1

        # Add mu1 center (we only vary mu2)
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'circle2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 1

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for center in centers:

            # Add mu2 center
            mu2['center']= 'np.array(['+str(center)+',0])'

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']


    # ==========================================================================
    #
    # Simulation 3: Squares moving closer
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two squares of 
    # equal diameter close to one another. For this reason, we vary the squares
    # center and, as usual, the number of subjects.
    #
    # ==========================================================================
    if simNo==3:

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # Add threshold c
        inputs['c'] = 2

        # Add mode
        inputs['mode'] = 1

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center (we only vary mu2)
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for center in centers:

            # Add mu2 center
            mu2['center']= 'np.array(['+str(center)+',0])'

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

    # ==========================================================================
    #
    # Simulation 4: Squares moving closer (but lower fwhm)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two squares of 
    # equal diameter close to one another. For this reason, we vary the squares
    # center and, as usual, the number of subjects. This differs from the 
    # previous setting in terms of smoothing and fwhm.
    #
    # ==========================================================================
    if simNo==4:

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # Add threshold c
        inputs['c'] = 2/3

        # Add mode
        inputs['mode'] = 1

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 1

        # Add mu1 center (we only vary mu2)
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 1

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for center in centers:

            # Add mu2 center
            mu2['center']= 'np.array(['+str(center)+',0])'

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

    # ==========================================================================
    #
    # Simulation 5: Squares moving closer (mode 2)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two squares of 
    # equal diameter close to one another. For this reason, we vary the squares
    # center and, as usual, the number of subjects.
    #
    # ==========================================================================
    if simNo==5:

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # Add threshold c
        inputs['c'] = 2

        # Add mode
        inputs['mode'] = 2

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center (we only vary mu2)
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for center in centers:

            # Add mu2 center
            mu2['center']= 'np.array(['+str(center)+',0])'

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

    # ==========================================================================
    #
    # Simulation 6: Squares moving closer (but lower fwhm) (mode 2)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two squares of 
    # equal diameter close to one another. For this reason, we vary the squares
    # center and, as usual, the number of subjects. This differs from the 
    # previous setting in terms of smoothing and fwhm.
    #
    # ==========================================================================
    if simNo==6:

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # Add threshold c
        inputs['c'] = 2/3

        # Add mode
        inputs['mode'] = 2

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 1

        # Add mu1 center (we only vary mu2)
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 1

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for center in centers:

            # Add mu2 center
            mu2['center']= 'np.array(['+str(center)+',0])'

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']


    # ==========================================================================
    #
    # Simulation 7: Squares moving closer (mode 3)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two squares of 
    # equal diameter close to one another. For this reason, we vary the squares
    # center and, as usual, the number of subjects.
    #
    # ==========================================================================
    if simNo==7:

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # Add threshold c
        inputs['c'] = 2

        # Add mode
        inputs['mode'] = 3

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center (we only vary mu2)
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for center in centers:

            # Add mu2 center
            mu2['center']= 'np.array(['+str(center)+',0])'

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

    # ==========================================================================
    #
    # Simulation 8: Squares moving closer (but lower fwhm) (mode 3)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two squares of 
    # equal diameter close to one another. For this reason, we vary the squares
    # center and, as usual, the number of subjects. This differs from the 
    # previous setting in terms of smoothing and fwhm.
    #
    # ==========================================================================
    if simNo==8:

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # Add threshold c
        inputs['c'] = 2/3

        # Add mode
        inputs['mode'] = 3

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 1

        # Add mu1 center (we only vary mu2)
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 1

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for center in centers:

            # Add mu2 center
            mu2['center']= 'np.array(['+str(center)+',0])'

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

    # ==========================================================================
    #
    # Simulation 9: Changing noise fwhm in one square (mode 3, high mu fwhm)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in changing the fwhm of the
    # noise in one of our fields but keeping the other constant (i.e. to see
    # if having different spatial correlations between the fields has an affect.
    # As always it is run across a range of subjects as well
    #
    # ==========================================================================
    if simNo==9:

        # These are our noise fwhms
        FWHM2s = np.linspace(1,6,26)

        # Add threshold c
        inputs['c'] = 2

        # Add mode
        inputs['mode'] = 3

        # Add FWHM for noise 1
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Save noise 1
        inputs['noise1'] = noise1

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3

        # Add mu2 center
        mu2['center'] = 'np.array([0,0])'

        # Add mu2 to inputs
        inputs['mu2'] = mu2

        # We will generate figures for these settings
        fg_FWHM2s = np.array([1,3,6])

        # Id for config file
        cfgId = 1

        # Loop through all FWHM2 settings
        for FWHM2 in FWHM2s:

            # Add FWHM for second noise field
            noise2['FWHM'] = '[0, ' + str(FWHM2) + ', ' + str(FWHM2) + ']'

            # Save noise field 2
            inputs['noise2'] = noise2   

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (FWHM2 in fg_FWHM2s):

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
        del inputs['noise2']['FWHM'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

    # ==========================================================================
    #
    # Simulation 10: Changing noise fwhm in one square (mode 3, low mu fwhm)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in changing the fwhm of the
    # noise in one of our fields but keeping the other constant (i.e. to see
    # if having different spatial correlations between the fields has an affect.
    # As always it is run across a range of subjects as well
    #
    # ==========================================================================
    if simNo==10:

        # These are our noise fwhms
        FWHM2s = np.linspace(1,6,26)

        # Add threshold c
        inputs['c'] = 2/3

        # Add mode
        inputs['mode'] = 3

        # Add FWHM for noise 1
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Save noises
        inputs['noise1'] = noise1

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 1

        # Add mu1 center
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 1

        # Add mu2 center
        mu2['center'] = 'np.array([0,0])'

        # Add mu2 to inputs
        inputs['mu2'] = mu2

        # We will generate figures for these settings
        fg_FWHM2s = np.array([1,3,6])

        # Id for config file
        cfgId = 1

        # Loop through all FWHM2 settings
        for FWHM2 in FWHM2s:

            # Add FWHM  for noise 2
            noise2['FWHM'] = '[0, ' + str(FWHM2) + ', ' + str(FWHM2) + ']'

            # Add noise 2 to inputs
            inputs['noise2'] = noise2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (FWHM2 in fg_FWHM2s):

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
        del inputs['noise2']['FWHM'], inputs['figGen'], inputs['cfgId'], inputs['nSub']


    # ==========================================================================
    #
    # Simulation 11: Squares moving closer (mode 3) (but heterogenous ramp on
    # one square)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two squares of 
    # equal diameter close to one another. For this reason, we vary the squares
    # center and, as usual, the number of subjects.
    #
    # ==========================================================================
    if simNo==11:

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # Add threshold c
        inputs['c'] = 2

        # Add mode
        inputs['mode'] = 3

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'heterogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center (we only vary mu2)
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for center in centers:

            # Add mu2 center
            mu2['center']= 'np.array(['+str(center)+',0])'

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

    # ==========================================================================
    #
    # Simulation 12: Squares moving closer (but lower fwhm) (mode 3)  (but 
    # heterogenous ramp on one square)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two squares of 
    # equal diameter close to one another. For this reason, we vary the squares
    # center and, as usual, the number of subjects. This differs from the 
    # previous setting in terms of smoothing and fwhm.
    #
    # ==========================================================================
    if simNo==12:

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # Add threshold c
        inputs['c'] = 2/3

        # Add mode
        inputs['mode'] = 3

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'heterogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 1

        # Add mu1 center (we only vary mu2)
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 1

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for center in centers:

            # Add mu2 center
            mu2['center']= 'np.array(['+str(center)+',0])'

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']


    # ==========================================================================
    #
    # Simulation 13: Squares moving closer (mode 3) (but heterogenous ramp on
    # both squares)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two squares of 
    # equal diameter close to one another. For this reason, we vary the squares
    # center and, as usual, the number of subjects.
    #
    # ==========================================================================
    if simNo==13:

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # Add threshold c
        inputs['c'] = 2

        # Add mode
        inputs['mode'] = 3

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'heterogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'heterogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center (we only vary mu2)
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for center in centers:

            # Add mu2 center
            mu2['center']= 'np.array(['+str(center)+',0])'

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

    # ==========================================================================
    #
    # Simulation 14: Squares moving closer (but lower fwhm) (mode 3)  (but 
    # heterogenous ramp on both square)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two squares of 
    # equal diameter close to one another. For this reason, we vary the squares
    # center and, as usual, the number of subjects. This differs from the 
    # previous setting in terms of smoothing and fwhm.
    #
    # ==========================================================================
    if simNo==14:

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # Add threshold c
        inputs['c'] = 2/3

        # Add mode
        inputs['mode'] = 3

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'heterogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'heterogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 1

        # Add mu1 center (we only vary mu2)
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 1

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for center in centers:

            # Add mu2 center
            mu2['center']= 'np.array(['+str(center)+',0])'

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']


    # ==========================================================================
    #
    # Simulation 15: Squares moving closer, one smaller than other
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two squares of 
    # equal diameter close to one another. For this reason, we vary the squares
    # center and, as usual, the number of subjects.
    #
    # ==========================================================================
    if simNo==15:

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # Add threshold c
        inputs['c'] = 2

        # Add mode
        inputs['mode'] = 3

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 20

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center (we only vary mu2)
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for center in centers:

            # Add mu2 center
            mu2['center']= 'np.array(['+str(center)+',0])'

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

    # ==========================================================================
    #
    # Simulation 16: Squares moving closer, one smaller than other (but lower fwhm)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two squares of 
    # equal diameter close to one another. For this reason, we vary the squares
    # center and, as usual, the number of subjects. This differs from the 
    # previous setting in terms of smoothing and fwhm.
    #
    # ==========================================================================
    if simNo==16:

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # Add threshold c
        inputs['c'] = 2/3

        # Add mode
        inputs['mode'] = 3

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 20

        # Add mu1 magnitude
        mu1['mag'] = 1

        # Add mu1 center (we only vary mu2)
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 1

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for center in centers:

            # Add mu2 center
            mu2['center']= 'np.array(['+str(center)+',0])'

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']


    # ==========================================================================
    #
    # Simulation 17: Varying correlation between fields
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we have two squares and we are varying the 
    # covariance between the noise in each field. 
    #
    # ==========================================================================
    if simNo==17:

        # These are our covariances
        noise_covs = np.arange(0.1,1.1,0.1)

        # Add threshold c
        inputs['c'] = 2

        # Add mode
        inputs['mode'] = 3

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 20

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3

        # Add mu2 center
        mu2['center'] = 'np.array([20,0])'

        # Add mu2 to inputs
        inputs['mu2'] = mu2

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for cov in noise_covs:

            # Add noise covariance
            inputs['noiseCov'] = cov

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['noiseCov'], inputs['figGen'], inputs['cfgId'], inputs['nSub']


    # ==========================================================================
    #
    # Simulation 18: Varying correlation between fields (lower SNR)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we have two squares and we are varying the 
    # covariance between the noise in each field. 
    #
    # ==========================================================================
    if simNo==18:

        # These are our covariances
        noise_covs = np.arange(0.1,1.1,0.1)

        # Add threshold c
        inputs['c'] = 2/3

        # Add mode
        inputs['mode'] = 3

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 20

        # Add mu1 magnitude
        mu1['mag'] = 1

        # Add mu1 center
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 1

        # Add mu2 center
        mu2['center'] = 'np.array([20,0])'

        # Add mu2 to inputs
        inputs['mu2'] = mu2

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for cov in noise_covs:

            # Add noise covariance
            inputs['noiseCov'] = cov

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['noiseCov'], inputs['figGen'], inputs['cfgId'], inputs['nSub']


    # ==========================================================================
    #
    # Simulation 19: Varying correlation between fields (higher SNR)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we have two squares and we are varying the 
    # covariance between the noise in each field. 
    #
    # ==========================================================================
    if simNo==19:

        # These are our covariances
        noise_covs = np.arange(0.1,1.1,0.1)

        # Add threshold c
        inputs['c'] = 2

        # Add mode
        inputs['mode'] = 3

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Save noises
        inputs['noise1'] = noise1
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 20

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center
        mu1['center'] = 'np.array([-20,0])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3

        # Add mu2 center
        mu2['center'] = 'np.array([20,0])'

        # Add mu2 to inputs
        inputs['mu2'] = mu2

        # We will generate figures for these settings
        fg_centers = np.array([-20,0,20,28])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for cov in noise_covs:

            # Add noise covariance
            inputs['noiseCov'] = cov

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and (center in fg_centers):

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
        del inputs['noiseCov'], inputs['figGen'], inputs['cfgId'], inputs['nSub']
    
    # Simulated over covariance range
    if simNo in [18,19]:

        # Variable to check if this is the first file weve looked at
        first = True

        # Loop through configuration files
        for cfgFile in cfgFiles:

            # ------------------------------------------------------------------
            # Load in inputs 
            # ------------------------------------------------------------------
            # Read in file
            with open(cfgFile, 'r') as stream:
                inputs = yaml.load(stream,Loader=yaml.FullLoader)

            # ------------------------------------------------------------------
            # Get directory of results for this config file
            # ------------------------------------------------------------------
            # Get configuration ID
            cfgId = int(inputs['cfgId'])

            # Results directory
            resDir = os.path.join(OutDir, 'sim'+str(simNo), 'cfg' + str(cfgId), 'RawResults')

            # ------------------------------------------------------------------
            # Get number of p values
            # ------------------------------------------------------------------ 
            # Get p values
            p = eval(inputs['p'])

            # Number of p values
            n_p = np.prod(p.shape)

            # ------------------------------------------------------------------
            # Currently there is a try except clause here. The reason for this
            # is that if a boundary is empty the code will error and no results
            # will be produced. In such a case, we skip that result.
            # ------------------------------------------------------------------ 
            try:

                # ------------------------------------------------------------------
                # Read in results
                # ------------------------------------------------------------------

                # Read in observed values for the estimated boundary. This will be a 
                # boolean array of ones and zeros representing observed violations
                # across simulations
                obs_est = pd.read_csv(os.path.join(resDir,'estSuccess.csv'), header=None, index_col=None)

                # Read in observed values for the true boundary. This will be a 
                # boolean array of ones and zeros representing observed violations
                # across simulations
                obs_true = pd.read_csv(os.path.join(resDir,'trueSuccess.csv'), header=None, index_col=None)

                # Read in observed values for the estimated boundary. This will be a 
                # boolean array of ones and zeros representing observed violations
                # across simulations (based on interpolation assessment)
                obs_est_intrp = pd.read_csv(os.path.join(resDir,'estSuccess_intrp.csv'), header=None, index_col=None)

                # Read in observed values for the true boundary. This will be a 
                # boolean array of ones and zeros representing observed violations
                # across simulations (based on interpolation assessment)
                obs_true_intrp = pd.read_csv(os.path.join(resDir,'trueSuccess_intrp.csv'), header=None, index_col=None)

                # ------------------------------------------------------------------
                # Get coverage probabilities
                # ------------------------------------------------------------------

                # Get the coverage probabilities from the observed results for the
                # estimated boundary
                covp_est = np.mean(obs_est.values,axis=0)[:]

                # Get the coverage probabilities from the observed results for the
                # true boundary
                covp_true = np.mean(obs_true.values,axis=0)[:]

                # Get the coverage probabilities from the observed results for the
                # estimated boundary (for coverage assessed using interpolation)
                covp_est_intrp = np.mean(obs_est_intrp.values,axis=0)[:]

                # Get the coverage probabilities from the observed results for the
                # true boundary (for coverage assessed using interpolation)
                covp_true_intrp = np.mean(obs_true_intrp.values,axis=0)[:]

                # ------------------------------------------------------------------
                # Get number of subjects and covariance between noise fields
                # ------------------------------------------------------------------

                # Number of subjects
                nSub = inputs['nSub']

                # Covariance between noise fields
                cov = inputs['noiseCov']

                # ------------------------------------------------------------------
                # Add coverage probabilities to table
                # ------------------------------------------------------------------
                # Line for table of estimated boundary results
                tableLine_est = np.concatenate((np.array([[cfgId,nSub,cov]]),\
                                                covp_est.reshape(1,n_p)),\
                                                axis=1)

                # Line for table of true boundary results
                tableLine_true = np.concatenate((np.array([[cfgId,nSub,cov]]),\
                                                 covp_true.reshape(1,n_p)),\
                                                 axis=1)

                # Line for table of estimated boundary interpolation assessed results
                tableLine_est_intrp = np.concatenate((np.array([[cfgId,nSub,cov]]),\
                                                      covp_est_intrp.reshape(1,n_p)),\
                                                      axis=1)

                # Line for table of true boundary interpolation assessed results
                tableLine_true_intrp = np.concatenate((np.array([[cfgId,nSub,cov]]),\
                                                       covp_true_intrp.reshape(1,n_p)),\
                                                       axis=1)

                # If this is the first cfg we've looked at, intialize the results tables
                if first:

                    # Initialize estimated boundary results table
                    table_est = pd.DataFrame(tableLine_est)

                    # Initialize true boundary results table
                    table_true = pd.DataFrame(tableLine_true)

                    # Initialize estimated boundary interpolated results table
                    table_est_intrp = pd.DataFrame(tableLine_est_intrp)

                    # Initialize true boundary interpolated results table
                    table_true_intrp = pd.DataFrame(tableLine_true_intrp)

                else:

                    # Append to existing estimated boundary results table
                    table_est = table_est.append(pd.DataFrame(tableLine_est))

                    # Append to existing true boundary results table
                    table_true = table_true.append(pd.DataFrame(tableLine_true))

                    # Append to existing estimated boundary interpolated results table
                    table_est_intrp = table_est_intrp.append(pd.DataFrame(tableLine_est_intrp))

                    # Append to existing true boundary interpolated results table
                    table_true_intrp = table_true_intrp.append(pd.DataFrame(tableLine_true_intrp))

                # ------------------------------------------------------------------
                # Get computation times
                # ------------------------------------------------------------------
                tableLine_time = pd.read_csv(os.path.join(resDir,'computationTime.csv'), header=None, index_col=None)

                # If this is the first cfg we've looked at, intialize the results tables
                if first:

                    # Initialize time table
                    table_time = pd.DataFrame(tableLine_time)

                else:

                    # Append to existing time table
                    table_time = table_time.append(pd.DataFrame(tableLine_time))

                # ------------------------------------------------------------------
                # Delete files
                # ------------------------------------------------------------------
                # Delete folder for this simulation
                shutil.rmtree(os.path.join(OutDir, 'sim'+str(simNo), 'cfg' + str(cfgId)))

                # We are no longer looking at the first configuration file
                if first:
                    first = False

            except:

                pass

        # ----------------------------------------------------------------------
        # Sort and save to csv
        # ----------------------------------------------------------------------
        # Make final results results directory
        fResDir = os.path.join(OutDir, 'sim'+str(simNo), 'FinalResults')
        if not os.path.exists(fResDir):
            os.mkdir(fResDir)

        # # Save times table
        # append_to_file(os.path.join(fResDir,'times.csv'), table_time)

        # Save estimated boundary results table
        append_to_file(os.path.join(fResDir,'estBdry.csv'), table_est)

        # Save true boundary results table
        append_to_file(os.path.join(fResDir,'trueBdry.csv'), table_true)

        # Save estimated boundary (with interpolation) results table
        append_to_file(os.path.join(fResDir,'estBdry_intrp.csv'), table_est_intrp)

        # Save true boundary (with interpolation) results table
        append_to_file(os.path.join(fResDir,'trueBdry_intrp.csv'), table_true_intrp)

        # ----------------------------------------------------------------------
        # Make figures
        # ----------------------------------------------------------------------

        # Column headers
        colhdr = ['cfgID', 'n', 'cov']+['p='+('%.2f' % p) for p in np.linspace(0,1,21)]

        # Assign column headers
        table_true_intrp.columns=colhdr
        table_est_intrp.columns=colhdr

        # List of n and p values
        n_values = np.unique(table_est_intrp['n'].values)
        c_values = np.unique(table_est_intrp['cov'].values)
        p_values = np.linspace(0,1,21)

        # Loop through all values of n
        for n in n_values:

            # Loop through all values of p
            for p in p_values:

                table_est_n = table_est_intrp[table_est_intrp['n']==n].sort_values('cov')
                table_true_n = table_true_intrp[table_true_intrp['n']==n].sort_values('cov')

                # Covariances
                cov_est_n = table_est_n[['cov']].values
                p_est_n = table_est_n[['p='+('%.2f' % p)]].values
                cov_true_n = table_true_n[['cov']].values
                p_true_n = table_true_n[['p='+('%.2f' % p)]].values

                plt.plot(cov_est_n,p_est_n,color="red",label="Estimated boundary")
                plt.plot(cov_true_n,p_true_n,color="blue",label="True boundary")
                plt.hlines(p, np.min(cov_est_n), np.max(cov_est_n),linestyles='dashed',label="Expected")

                # Title
                plt.title("Coverage (" + str(np.int(100*p)) + "% probability, " + str(int(n)) + " subjects)")

                # Axes
                plt.xlabel("Covariance between noise fields")
                plt.ylabel("Observed coverage")

                # Make axis a bit clearer
                plt.ylim((np.min(p_true_n)-0.02,1))
                
                # Legend
                plt.legend()

                # Save plots
                plt.savefig(os.path.join(fResDir, 'cov_vs_obsp_truep'+str(np.int(100*p))+'_n'+str(np.int(n))+'.png'))

                # Clear figure
                plt.clf()

        # Loop through all values of fwhm2
        for c in c_values:

            # Loop through all values of p
            for p in p_values:

                table_est_c = table_est_intrp[table_est_intrp['cov']==c].sort_values('n')
                table_true_c = table_true_intrp[table_true_intrp['cov']==c].sort_values('n')

                # n and p for this covariance
                n_est_c = table_est_c[['n']].values
                p_est_c = table_est_c[['p='+('%.2f' % p)]].values
                n_true_c = table_true_c[['n']].values
                p_true_c = table_true_c[['p='+('%.2f' % p)]].values

                plt.plot(n_est_c,p_est_c,color="red",label="Estimated boundary")
                plt.plot(n_true_c,p_true_c,color="blue",label="True boundary")
                plt.hlines(p, np.min(n_est_c), np.max(n_est_c),linestyles='dashed',label="Expected")

                # Title
                plt.title("Coverage (" + str(np.int(100*p)) + "% probability, noise covariance " + ('%.2f' % c) + ")")

                # Axes
                plt.xlabel("Number of subjects")
                plt.ylabel("Observed coverage")
                
                # Make axis a bit clearer
                plt.ylim((np.min(p_true_c)-0.02,1))
                
                # Legend
                plt.legend()

                # Save plots
                plt.savefig(os.path.join(fResDir, 'n_vs_obsp_truep'+str(np.int(100*p))+'_cov'+('%.2f' % c)+'.png'))

                # Clear figure
                plt.clf()
    #--------------------------------------------------------------------------------------
    # Save the baseline configuration (Note: This must be the last file output as the bash
    # script used for running simulations on the cluster will take this files existence as
    # a sign to run the next stage of the simulations)
    #--------------------------------------------------------------------------------------
    # Save the yml
    with open(os.path.join(simDir,'cfgs','baseline_cfg.yml'), 'w') as outfile:
        yaml.dump(inputs, outfile, default_flow_style=False)



#generateCfgs('/home/tommaullin/Documents/ConfRes/tmp/sim12', 12)
