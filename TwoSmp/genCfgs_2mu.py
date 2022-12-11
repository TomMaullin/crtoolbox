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
    inputs['nReals'] = 2500

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
        inputs['c'] = 1/2

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
        mu1['type'] = 'circle2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3/4

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
        mu2['mag'] = 3/4

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
        inputs['c'] = 1/2

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
        mu1['mag'] = 3/4

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
        mu2['mag'] = 3/4

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
        inputs['c'] = 1/2

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
        mu1['mag'] = 3/4

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
        mu2['mag'] = 3/4

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
        inputs['c'] = 1/2

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
        mu1['mag'] = 3/4

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
        mu2['mag'] = 3/4

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
        inputs['c'] = 1/2

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
        mu1['mag'] = 3/4

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
        mu2['mag'] = 3/4

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
        inputs['c'] = 1/2

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
        mu1['mag'] = 3/4

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
        mu2['mag'] = 3/4

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
        inputs['c'] = 1/2

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
        mu1['mag'] = 3/4

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
        mu2['mag'] = 3/4

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
        inputs['c'] = 1/2

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
        mu1['mag'] = 3/4

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
        mu2['mag'] = 3/4

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
    # Simulation 17: Varying correlation between fields (higher SNR)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we have two squares and we are varying the 
    # covariance between the noise in each field. 
    #
    # ==========================================================================
    if simNo==17:

        # These are our covariances
        noise_corrs = np.arange(-1.0,1.1,0.1)

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
        fg_corrs = np.array([-0.8,0,0.8])

        # Id for config file
        cfgId = 1

        # Loop through all correlation settings
        for corr in noise_corrs:

            # Add noise correlation
            inputs['noiseCorr'] = str(corr)

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and np.any(np.isclose(fg_corrs,corr)):

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
        del inputs['noiseCorr'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

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
        noise_corrs = np.arange(-1.0,1.1,0.1)

        # Add threshold c
        inputs['c'] = 1/2

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
        mu1['mag'] = 3/4

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
        mu2['mag'] = 3/4

        # Add mu2 center
        mu2['center'] = 'np.array([20,0])'

        # Add mu2 to inputs
        inputs['mu2'] = mu2

        # We will generate figures for these settings
        fg_corrs = np.array([-0.8,0,0.8])

        # Id for config file
        cfgId = 1

        # Loop through all correlation settings
        for corr in noise_corrs:

            # Add noise correlation
            inputs['noiseCorr'] = str(corr)

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and np.any(np.isclose(fg_corrs,corr)):

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
        del inputs['noiseCorr'], inputs['figGen'], inputs['cfgId'], inputs['nSub']



    # ==========================================================================
    #
    # Simulation 19: Ramps moving closer 
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two ramps close
    # to one another. We vary the ramps slopes and, as usual, the number of
    # subjects.
    #
    # ==========================================================================
    if simNo==19:

        # These are our ramp gradients
        grads = 4*np.arange(0.5,3.6,0.1)

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


        # We will generate figures for these settings
        fg_grads = 4*np.array([1,2,3])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for grad in grads:

            # Add mu1 type
            mu1['type'] = 'ramp2D'

            # Add mu1 fwhm
            mu1['orient'] = 'horizontal'

            # Add mu1 a
            mu1['a'] = '%.15f' % (2-grad)

            # Add mu1 b
            mu1['b'] = '%.15f' % (2+grad)

            # Add mu1 to inputs
            inputs['mu1'] = mu1

            # Add mu2 type
            mu2['type'] = 'ramp2D'

            # Add mu2 fwhm
            mu2['orient'] = 'vertical'

            # Add mu2 a
            mu2['a'] = '%.15f' % (2-grad)

            # Add mu2 b
            mu2['b'] = '%.15f' % (2+grad)

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                print(grad, fg_grads, grad in fg_grads)
                if (nSub in fg_nSubs) and np.any(np.isclose(fg_grads,grad)):

                    # In this case we do want to save  figures
                    inputs['figGen']=1

                else:

                    # In this case we do want to save  figures
                    inputs['figGen']=0

                print(inputs['figGen'])

                # Save the yml
                with open(os.path.join(simDir,'cfgs','cfg'+str(cfgId)+'.yml'), 'w') as outfile:
                    yaml.dump(inputs, outfile, default_flow_style=False)

                # Incremement cfgID
                cfgId = cfgId + 1

        # Delete fields which vary across simulation
        del inputs['mu1']['a'], inputs['mu1']['b'], inputs['mu2']['a'], inputs['mu2']['b'], inputs['figGen'], inputs['cfgId'], inputs['nSub']


    # ==========================================================================
    #
    # Simulation 20: Ramps moving closer (but lower fwhm)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two ramps close
    # to one another. We vary the ramps slopes and, as usual, the number of
    # subjects. This differs from the previous setting in terms of smoothing
    # and fwhm.
    #
    # ==========================================================================
    if simNo==20:

        # These are our ramp gradients
        grads = np.arange(0.5,3.6,0.1)

        # Add threshold c
        inputs['c'] = 1/2

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


        # We will generate figures for these settings
        fg_grads = np.array([1,2,3])

        # Id for config file
        cfgId = 1

        # Loop through all center settings
        for grad in grads:

            # Add mu1 type
            mu1['type'] = 'ramp2D'

            # Add mu1 fwhm
            mu1['orient'] = 'horizontal'

            # Add mu1 a
            mu1['a'] = '%.15f' % (1/2-grad)

            # Add mu1 b
            mu1['b'] = '%.15f' % (1/2+grad)

            # Add mu1 to inputs
            inputs['mu1'] = mu1

            # Add mu2 type
            mu2['type'] = 'ramp2D'

            # Add mu2 fwhm
            mu2['orient'] = 'vertical'

            # Add mu2 a
            mu2['a'] = '%.15f' % (1/2-grad)

            # Add mu2 b
            mu2['b'] = '%.15f' % (1/2+grad)

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and np.any(np.isclose(fg_grads,grad)):

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
        del inputs['mu1']['a'], inputs['mu1']['b'], inputs['mu2']['a'], inputs['mu2']['b'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

    # ==========================================================================
    #
    # Simulation 21: Increasing noise magnitude
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in increasing the magnitude
    # of one of the noise fields. We also vary the number of subjects.
    #
    # ==========================================================================
    if simNo==21:

        # These are our noise 2 magnitudes
        mags = np.arange(1,3.2,0.2)

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

        # Add mu2 center
        mu2['center']= 'np.array([20,0])'

        # Add mu2 to inputs
        inputs['mu2'] = mu2

        # We will generate figures for these settings
        fg_mags = np.array([1,1.5,2])

        # Id for config file
        cfgId = 1

        # Loop through all noise magnitude settings
        for mag in mags:

            # Change noise2 mag
            noise2['mag']=str(mag)

            # Add noise2 to inputs
            inputs['noise2'] = noise2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and np.any(np.isclose(fg_mags,mag)):

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
        del inputs['noise2']['mag'], inputs['figGen'], inputs['cfgId'], inputs['nSub']
        
    # ==========================================================================
    #
    # Simulation 22: Increasing noise magnitude (but lower SNR)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in increasing the magnitude
    # of one of the noise fields. We also vary the number of subjects.This
    # differs from the previous setting in terms of smoothing and fwhm.
    #
    # ==========================================================================
    if simNo==22:

        # These are our noise 2 magnitudes
        mags = np.arange(1,3.2,0.2)

        # Add threshold c
        inputs['c'] = 1/2

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

        # Save noise 1
        inputs['noise1'] = noise1

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3/4

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
        mu2['mag'] = 3/4

        # Add mu2 center
        mu2['center']= 'np.array([20,0])'

        # Add mu2 to inputs
        inputs['mu2'] = mu2

        # We will generate figures for these settings
        fg_mags = np.array([1,1.5,2])

        # Id for config file
        cfgId = 1

        # Loop through all noise magnitude settings
        for mag in mags:

            # Change noise2 mag
            noise2['mag']=str(mag)

            # Add noise2 to inputs
            inputs['noise2'] = noise2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and np.any(np.isclose(fg_mags,mag)):

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
        del inputs['noise2']['mag'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

    # --------------------------------------------------------------------------
    # Reviewer Requested Simulations
    # --------------------------------------------------------------------------
    # The simulations that follow were added to this file during the first week
    # of December 2022 following reviewer requested feedback.
    # --------------------------------------------------------------------------

    # ==========================================================================
    #
    # Simulation 23: Varying smoothness of signal
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in varying the smoothess of 
    # the signal, with particular interest being paid to the setting of zero
    # spatial correlation.
    #
    # ==========================================================================
    if simNo==23:

        # These are our signal smoothness values
        smooths = np.arange(0,8,0.2)

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

        # Save noise 1
        inputs['noise1'] = noise1

        # Save noise 2
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center (we only vary mu2)
        mu1['center'] = 'np.array([-20,0])'

        # Add mu2 type
        mu2['type'] = 'square2D' 


        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3

        # Add mu2 center
        mu2['center']= 'np.array([20,0])'

        # We will generate figures for these settings
        fg_smooths = np.array([0,2,4,6])

        # Id for config file
        cfgId = 1

        # Loop through all noise magnitude settings
        for smooth in smooths:

            # Add mu1 fwhm
            mu1['fwhm'] = 'np.array([' + str(smooth) + ',' + str(smooth) +'])'

            # Add mu2 fwhm
            mu2['fwhm'] = 'np.array([' + str(smooth) + ',' + str(smooth) +'])'

            # Add mu1 to inputs
            inputs['mu1'] = mu1

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and np.any(np.isclose(fg_smooths,smooth)):

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
        del inputs['mu1']['fwhm'], inputs['mu2']['fwhm'], inputs['figGen'], inputs['cfgId'], inputs['nSub']
   

    # ==========================================================================
    #
    # Simulation 24: Varying smoothness of signal (lower SNR)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in varying the smoothess of 
    # the signal, with particular interest being paid to the setting of zero
    # spatial correlation.
    #
    # ==========================================================================
    if simNo==24:

        # These are our signal smoothness values
        smooths = np.arange(0,8,0.2)

        # Add threshold c
        inputs['c'] = 1/2

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

        # Save noise 1
        inputs['noise1'] = noise1

        # Save noise 2
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3/4

        # Add mu1 center
        mu1['center'] = 'np.array([-20,0])'

        # Add mu2 type
        mu2['type'] = 'square2D' 


        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3/4

        # Add mu2 center
        mu2['center']= 'np.array([20,0])'

        # We will generate figures for these settings
        fg_smooths = np.array([0,2,4,6])

        # Id for config file
        cfgId = 1

        # Loop through all noise magnitude settings
        for smooth in smooths:

            # Add mu1 fwhm
            mu1['fwhm'] = 'np.array([' + str(smooth) + ',' + str(smooth) +'])'

            # Add mu2 fwhm
            mu2['fwhm'] = 'np.array([' + str(smooth) + ',' + str(smooth) +'])'

            # Add mu1 to inputs
            inputs['mu1'] = mu1

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and np.any(np.isclose(fg_smooths,smooth)):

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
        del inputs['mu1']['fwhm'], inputs['mu2']['fwhm'], inputs['figGen'], inputs['cfgId'], inputs['nSub']


    # ==========================================================================
    #
    # Simulation 25: Varying smoothness of noise 
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in varying the smoothess of 
    # the noise, with particular interest being paid to the setting of zero
    # spatial correlation.
    #
    # ==========================================================================
    if simNo==25:

        # These are our signal smoothness values
        smooths = np.arange(0,8,0.2)

        # Add threshold c
        inputs['c'] = 2

        # Add mode
        inputs['mode'] = 3

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center
        mu1['center'] = 'np.array([-20,0])'

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3

        # Add mu2 center
        mu2['center']= 'np.array([20,0])'

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 to inputs
        inputs['mu2'] = mu2

        # We will generate figures for these settings
        fg_smooths = np.array([0,2,4,6])

        # Id for config file
        cfgId = 1

        # Loop through all noise magnitude settings
        for smooth in smooths:

            # Add FWHM for noise
            noise1['FWHM'] = '[0, ' + str(smooth) + ', ' + str(smooth) + ']'

            # Add FWHM  for noise 2
            noise2['FWHM'] = '[0, ' + str(smooth) + ', ' + str(smooth) + ']'

            # Save noise 1
            inputs['noise1'] = noise1

            # Save noise 2
            inputs['noise2'] = noise2


            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and np.any(np.isclose(fg_smooths,smooth)):

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
        del inputs['noise1']['FWHM'], inputs['noise2']['FWHM'], inputs['figGen'], inputs['cfgId'], inputs['nSub']
   

    # ==========================================================================
    #
    # Simulation 26: Varying smoothness of noise (lower SNR)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in varying the smoothess of 
    # the noise, with particular interest being paid to the setting of zero
    # spatial correlation.
    #
    # ==========================================================================
    if simNo==26:

        # These are our signal smoothness values
        smooths = np.arange(0,8,0.2)

        # Add threshold c
        inputs['c'] = 1/2

        # Add mode
        inputs['mode'] = 3

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3/4

        # Add mu1 center
        mu1['center'] = 'np.array([-20,0])'

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3/4

        # Add mu2 center
        mu2['center']= 'np.array([20,0])'

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 to inputs
        inputs['mu2'] = mu2

        # We will generate figures for these settings
        fg_smooths = np.array([0,2,4,6])

        # Id for config file
        cfgId = 1

        # Loop through all noise magnitude settings
        for smooth in smooths:

            # Add FWHM for noise
            noise1['FWHM'] = '[0, ' + str(smooth) + ', ' + str(smooth) + ']'

            # Add FWHM  for noise 2
            noise2['FWHM'] = '[0, ' + str(smooth) + ', ' + str(smooth) + ']'

            # Save noise 1
            inputs['noise1'] = noise1

            # Save noise 2
            inputs['noise2'] = noise2


            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and np.any(np.isclose(fg_smooths,smooth)):

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
        del inputs['noise1']['FWHM'], inputs['noise2']['FWHM'], inputs['figGen'], inputs['cfgId'], inputs['nSub']



    # ==========================================================================
    #
    # Simulation 27: Varying c
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in varying the threhold c.
    #
    # ==========================================================================
    if simNo==27:

        # These are our signal smoothness values
        cs = np.arange(0,4,0.2)

        # Add mode
        inputs['mode'] = 3

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Save noise 1
        inputs['noise1'] = noise1

        # Save noise 2
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center
        mu1['center'] = 'np.array([-20,0])'

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3

        # Add mu2 center
        mu2['center']= 'np.array([20,0])'

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 to inputs
        inputs['mu2'] = mu2

        # We will generate figures for these settings
        fg_cs = np.array([0,1,2,3])

        # Id for config file
        cfgId = 1

        # Loop through all noise magnitude settings
        for c in cs:

            # Add threshold c
            inputs['c'] = str(c)

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and np.any(np.isclose(fg_cs,c)):

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
        del inputs['c'], inputs['figGen'], inputs['cfgId'], inputs['nSub']



    # ==========================================================================
    #
    # Simulation 28: Varying c (Low SNR)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in varying the threhold c.
    #
    # ==========================================================================
    if simNo==28:

        # These are our signal smoothness values
        cs = np.arange(0,4,0.2)/4

        # Add mode
        inputs['mode'] = 3

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Save noise 1
        inputs['noise1'] = noise1

        # Save noise 2
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'square2D' 

        # Add mu1 radius
        mu1['r'] = 30

        # Add mu1 magnitude
        mu1['mag'] = 3/4

        # Add mu1 center
        mu1['center'] = 'np.array([-20,0])'

        # Add mu2 type
        mu2['type'] = 'square2D' 

        # Add mu2 radius
        mu2['r'] = 30

        # Add mu2 magnitude
        mu2['mag'] = 3/4

        # Add mu2 center
        mu2['center']= 'np.array([20,0])'

        # Add mu1 fwhm
        mu1['fwhm'] = 'np.array([5,5])'

        # Add mu2 fwhm
        mu2['fwhm'] = 'np.array([5,5])'

        # Add mu1 to inputs
        inputs['mu1'] = mu1

        # Add mu2 to inputs
        inputs['mu2'] = mu2

        # We will generate figures for these settings
        fg_cs = np.array([0,1,2,3])/4

        # Id for config file
        cfgId = 1

        # Loop through all noise magnitude settings
        for c in cs:

            # Add threshold c
            inputs['c'] = str(c)

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and np.any(np.isclose(fg_cs,c)):

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
        del inputs['c'], inputs['figGen'], inputs['cfgId'], inputs['nSub']



    # ==========================================================================
    #
    # Simulation 29: 3 circles, intersection of CRs
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in comparing the method to
    # just intersecting CRs. A specially generated signal of three circles 
    # diagonal is used (one for each condition and one shared) with varying 
    # signal smoothness, and we run the intersection of CRs method here.
    #
    # ==========================================================================
    if simNo==29:

        # These are our signal smoothness values
        smooths = np.arange(0,8,0.2)

        # Create muBoth specification
        muBoth = {}

        # Add mode
        inputs['mode'] = 3

        inputs['Seperate'] = 1

        # Add threshold c
        inputs['c'] = 2

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Save noise 1
        inputs['noise1'] = noise1

        # Save noise 2
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'circle2D' 

        # Add mu1 radius
        mu1['r'] = 12

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center
        mu1['center'] = 'np.array([-30,-30])'

        # Add mu2 type
        mu2['type'] = 'circle2D' 

        # Add mu2 radius
        mu2['r'] = 12

        # Add mu2 magnitude
        mu2['mag'] = 3

        # Add mu2 center
        mu2['center']= 'np.array([30,30])'

        # Add muBoth type
        muBoth['type'] = 'circle2D' 

        # Add muBoth radius
        muBoth['r'] = 12

        # Add muBoth magnitude
        muBoth['mag'] = 3

        # Add muBoth center
        muBoth['center']= 'np.array([0,0])'

        # We will generate figures for these settings
        fg_smooths = np.array([0,2,4,6])

        # Id for config file
        cfgId = 1

        # Loop through all noise magnitude settings
        for smooth in smooths:

            # Add mu1 fwhm
            mu1['fwhm'] = 'np.array([' + str(2+smooth) + ',' + str(2+smooth) +'])'

            # Add mu2 fwhm
            mu2['fwhm'] = 'np.array([' + str(2+smooth) + ',' + str(2+smooth) +'])'

            # Add muBoth fwhm
            muBoth['fwhm'] = 'np.array([' + str(10-smooth) + ',' + str(10-smooth) +'])'

            # Add mu1 to inputs
            inputs['mu1'] = mu1

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Add muBoth to inputs
            inputs['muBoth'] = muBoth

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and np.any(np.isclose(fg_smooths,smooth)):

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
        del inputs['mu1']['fwhm'], inputs['mu2']['fwhm'], inputs['muBoth']['fwhm'], inputs['figGen'], inputs['cfgId'], inputs['nSub']


    # ==========================================================================
    #
    # Simulation 30: 3 circles, intersection of CRs (Low SNR)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in comparing the method to
    # just intersecting CRs. A specially generated signal of three circles 
    # diagonal is used (one for each condition and one shared) with varying 
    # signal smoothness, and we run the intersection of CRs method here.
    #
    # ==========================================================================
    if simNo==30:

        # These are our signal smoothness values
        smooths = np.arange(0,8,0.2)

        # Create muBoth specification
        muBoth = {}

        # Add mode
        inputs['mode'] = 3

        inputs['Seperate'] = 1

        # Add threshold c
        inputs['c'] = 1/2

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Save noise 1
        inputs['noise1'] = noise1

        # Save noise 2
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'circle2D' 

        # Add mu1 radius
        mu1['r'] = 12

        # Add mu1 magnitude
        mu1['mag'] = 3/4

        # Add mu1 center
        mu1['center'] = 'np.array([-30,-30])'

        # Add mu2 type
        mu2['type'] = 'circle2D' 

        # Add mu2 radius
        mu2['r'] = 12

        # Add mu2 magnitude
        mu2['mag'] = 3/4

        # Add mu2 center
        mu2['center']= 'np.array([30,30])'

        # Add muBoth type
        muBoth['type'] = 'circle2D' 

        # Add muBoth radius
        muBoth['r'] = 12

        # Add muBoth magnitude
        muBoth['mag'] = 3/4

        # Add muBoth center
        muBoth['center']= 'np.array([0,0])'

        # We will generate figures for these settings
        fg_smooths = np.array([0,2,4,6])

        # Id for config file
        cfgId = 1

        # Loop through all noise magnitude settings
        for smooth in smooths:

            # Add mu1 fwhm
            mu1['fwhm'] = 'np.array([' + str(2+smooth) + ',' + str(2+smooth) +'])'

            # Add mu2 fwhm
            mu2['fwhm'] = 'np.array([' + str(2+smooth) + ',' + str(2+smooth) +'])'

            # Add muBoth fwhm
            muBoth['fwhm'] = 'np.array([' + str(10-smooth) + ',' + str(10-smooth) +'])'

            # Add mu1 to inputs
            inputs['mu1'] = mu1

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Add muBoth to inputs
            inputs['muBoth'] = muBoth

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and np.any(np.isclose(fg_smooths,smooth)):

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
        del inputs['mu1']['fwhm'], inputs['mu2']['fwhm'], inputs['muBoth']['fwhm'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

    # ==========================================================================
    #
    # Simulation 31: 3 circles, intersection of CRs (conjunction method)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in comparing the method to
    # just intersecting CRs. A specially generated signal of three circles 
    # diagonal is used (one for each condition and one shared) with varying 
    # signal smoothness, and we run the conjunction method here.
    #
    # ==========================================================================
    if simNo==31:

        # These are our signal smoothness values
        smooths = np.arange(0,8,0.2)

        # Create muBoth specification
        muBoth = {}

        # Add mode
        inputs['mode'] = 3

        # Add threshold c
        inputs['c'] = 2

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Save noise 1
        inputs['noise1'] = noise1

        # Save noise 2
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'circle2D' 

        # Add mu1 radius
        mu1['r'] = 12

        # Add mu1 magnitude
        mu1['mag'] = 3

        # Add mu1 center
        mu1['center'] = 'np.array([-30,-30])'

        # Add mu2 type
        mu2['type'] = 'circle2D' 

        # Add mu2 radius
        mu2['r'] = 12

        # Add mu2 magnitude
        mu2['mag'] = 3

        # Add mu2 center
        mu2['center']= 'np.array([30,30])'

        # Add muBoth type
        muBoth['type'] = 'circle2D' 

        # Add muBoth radius
        muBoth['r'] = 12

        # Add muBoth magnitude
        muBoth['mag'] = 3

        # Add muBoth center
        muBoth['center']= 'np.array([0,0])'

        # We will generate figures for these settings
        fg_smooths = np.array([0,2,4,6])

        # Id for config file
        cfgId = 1

        # Loop through all noise magnitude settings
        for smooth in smooths:

            # Add mu1 fwhm
            mu1['fwhm'] = 'np.array([' + str(2+smooth) + ',' + str(2+smooth) +'])'

            # Add mu2 fwhm
            mu2['fwhm'] = 'np.array([' + str(2+smooth) + ',' + str(2+smooth) +'])'

            # Add muBoth fwhm
            muBoth['fwhm'] = 'np.array([' + str(10-smooth) + ',' + str(10-smooth) +'])'

            # Add mu1 to inputs
            inputs['mu1'] = mu1

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Add muBoth to inputs
            inputs['muBoth'] = muBoth

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and np.any(np.isclose(fg_smooths,smooth)):

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
        del inputs['mu1']['fwhm'], inputs['mu2']['fwhm'], inputs['muBoth']['fwhm'], inputs['figGen'], inputs['cfgId'], inputs['nSub']


    # ==========================================================================
    #
    # Simulation 32: 3 circles, intersection of CRs (conjunction method) (Low 
    # SNR)
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in comparing the method to
    # just intersecting CRs. A specially generated signal of three circles 
    # diagonal is used (one for each condition and one shared) with varying 
    # signal smoothness, and we run the conjunction method here.
    #
    # ==========================================================================
    if simNo==32:

        # These are our signal smoothness values
        smooths = np.arange(0,8,0.2)

        # Create muBoth specification
        muBoth = {}

        # Add mode
        inputs['mode'] = 3

        # Add threshold c
        inputs['c'] = 1/2

        # Add type for noise 1
        noise1['type'] = 'homogen'

        # Add type for noise 2
        noise2['type'] = 'homogen'

        # Add FWHM for noise
        noise1['FWHM'] = '[0, 3, 3]'

        # Add FWHM  for noise 2
        noise2['FWHM'] = '[0, 3, 3]'

        # Save noise 1
        inputs['noise1'] = noise1

        # Save noise 2
        inputs['noise2'] = noise2

        # Add mu1 type
        mu1['type'] = 'circle2D' 

        # Add mu1 radius
        mu1['r'] = 12

        # Add mu1 magnitude
        mu1['mag'] = 3/4

        # Add mu1 center
        mu1['center'] = 'np.array([-30,-30])'

        # Add mu2 type
        mu2['type'] = 'circle2D' 

        # Add mu2 radius
        mu2['r'] = 12

        # Add mu2 magnitude
        mu2['mag'] = 3/4

        # Add mu2 center
        mu2['center']= 'np.array([30,30])'

        # Add muBoth type
        muBoth['type'] = 'circle2D' 

        # Add muBoth radius
        muBoth['r'] = 12

        # Add muBoth magnitude
        muBoth['mag'] = 3/4

        # Add muBoth center
        muBoth['center']= 'np.array([0,0])'

        # We will generate figures for these settings
        fg_smooths = np.array([0,2,4,6])

        # Id for config file
        cfgId = 1

        # Loop through all noise magnitude settings
        for smooth in smooths:

            # Add mu1 fwhm
            mu1['fwhm'] = 'np.array([' + str(2+smooth) + ',' + str(2+smooth) +'])'

            # Add mu2 fwhm
            mu2['fwhm'] = 'np.array([' + str(2+smooth) + ',' + str(2+smooth) +'])'

            # Add muBoth fwhm
            muBoth['fwhm'] = 'np.array([' + str(10-smooth) + ',' + str(10-smooth) +'])'

            # Add mu1 to inputs
            inputs['mu1'] = mu1

            # Add mu2 to inputs
            inputs['mu2'] = mu2

            # Add muBoth to inputs
            inputs['muBoth'] = muBoth

            # Loop through all nSub settings
            for nSub in nSubs:

                # Add nSub to inputs
                inputs['nSub'] = int(nSub)

                # Save cfg ID (handy to have around)
                inputs['cfgId'] = int(cfgId)

                # Record if we want to save figures for this design or not
                if (nSub in fg_nSubs) and np.any(np.isclose(fg_smooths,smooth)):

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
        del inputs['mu1']['fwhm'], inputs['mu2']['fwhm'], inputs['muBoth']['fwhm'], inputs['figGen'], inputs['cfgId'], inputs['nSub']



    

    #--------------------------------------------------------------------------------------
    # Save the baseline configuration (Note: This must be the last file output as the bash
    # script used for running simulations on the cluster will take this files existence as
    # a sign to run the next stage of the simulations)
    #--------------------------------------------------------------------------------------
    # Save the yml
    with open(os.path.join(simDir,'cfgs','baseline_cfg.yml'), 'w') as outfile:
        yaml.dump(inputs, outfile, default_flow_style=False)

#generateCfgs('/home/tommaullin/Documents/ConfRes/tmp/sim19', 19)
