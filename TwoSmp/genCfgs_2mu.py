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

        # These are our sample sizes:
        nSubs = np.linspace(40,500,24)

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # New empty inputs directory
        inputs={}

        # Add output directory
        inputs['OutDir'] = OutDir

        # Add number of realizations
        inputs['nReals'] = 1000

        # Add threshold c
        inputs['c'] = 2

        # Add mode
        inputs['mode'] = 1

        # Add range of p-values
        inputs['p'] = 'np.linspace(0,1,21)'

        # Add simulation number (helps identify simulation)
        inputs['simNo'] = simNo

        # Add FWHM
        inputs['FWHM'] = '[0, 3, 3]'

        # Add number of bootstraps
        inputs['nBoot'] = 10000

        # Add tau expression
        inputs['tau'] = '1/np.sqrt(nSub)'

        # Create mu1 specification
        mu1 = {}

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

        # Create mu2 specification
        mu2 = {}

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
        fg_nSubs = np.array([100,300,500])

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

        #--------------------------------------------------------------------------------------
        # Save the baseline configuration (Note: This must be the last file output as the bash
        # script used for running simulations on the cluster will take this files existence as
        # a sign to run the next stage of the simulations)
        #--------------------------------------------------------------------------------------
        # Delete fields which vary across simulation
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

        # Save the yml
        with open(os.path.join(simDir,'cfgs','baseline_cfg.yml'), 'w') as outfile:
            yaml.dump(inputs, outfile, default_flow_style=False)

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

        # These are our sample sizes:
        nSubs = np.linspace(40,500,24)

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # New empty inputs directory
        inputs={}

        # Add output directory
        inputs['OutDir'] = OutDir

        # Add number of realizations
        inputs['nReals'] = 500

        # Add threshold c
        inputs['c'] = 2/3

        # Add mode
        inputs['mode'] = 1

        # Add range of p-values
        inputs['p'] = 'np.linspace(0,1,21)'

        # Add simulation number (helps identify simulation)
        inputs['simNo'] = simNo

        # Add FWHM
        inputs['FWHM'] = '[0, 3, 3]'

        # Add number of bootstraps
        inputs['nBoot'] = 5000

        # Add tau expression
        inputs['tau'] = '1/np.sqrt(nSub)'

        # Create mu1 specification
        mu1 = {}

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

        # Create mu2 specification
        mu2 = {}

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
        fg_nSubs = np.array([100,300,500])

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

        #--------------------------------------------------------------------------------------
        # Save the baseline configuration (Note: This must be the last file output as the bash
        # script used for running simulations on the cluster will take this files existence as
        # a sign to run the next stage of the simulations)
        #--------------------------------------------------------------------------------------
        # Delete fields which vary across simulation
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

        # Save the yml
        with open(os.path.join(simDir,'cfgs','baseline_cfg.yml'), 'w') as outfile:
            yaml.dump(inputs, outfile, default_flow_style=False)


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

        # These are our sample sizes:
        nSubs = np.linspace(40,500,24)

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # New empty inputs directory
        inputs={}

        # Add output directory
        inputs['OutDir'] = OutDir

        # Add number of realizations
        inputs['nReals'] = 1000

        # Add threshold c
        inputs['c'] = 2

        # Add mode
        inputs['mode'] = 1

        # Add range of p-values
        inputs['p'] = 'np.linspace(0,1,21)'

        # Add simulation number (helps identify simulation)
        inputs['simNo'] = simNo

        # Add FWHM
        inputs['FWHM'] = '[0, 3, 3]'

        # Add number of bootstraps
        inputs['nBoot'] = 10000

        # Add tau expression
        inputs['tau'] = '1/np.sqrt(nSub)'

        # Create mu1 specification
        mu1 = {}

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

        # Create mu2 specification
        mu2 = {}

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
        fg_nSubs = np.array([100,300,500])

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

        #--------------------------------------------------------------------------------------
        # Save the baseline configuration (Note: This must be the last file output as the bash
        # script used for running simulations on the cluster will take this files existence as
        # a sign to run the next stage of the simulations)
        #--------------------------------------------------------------------------------------
        # Delete fields which vary across simulation
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

        # Save the yml
        with open(os.path.join(simDir,'cfgs','baseline_cfg.yml'), 'w') as outfile:
            yaml.dump(inputs, outfile, default_flow_style=False)

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

        # These are our sample sizes:
        nSubs = np.linspace(40,500,24)

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # New empty inputs directory
        inputs={}

        # Add output directory
        inputs['OutDir'] = OutDir

        # Add number of realizations
        inputs['nReals'] = 500

        # Add threshold c
        inputs['c'] = 2/3

        # Add mode
        inputs['mode'] = 1

        # Add range of p-values
        inputs['p'] = 'np.linspace(0,1,21)'

        # Add simulation number (helps identify simulation)
        inputs['simNo'] = simNo

        # Add FWHM
        inputs['FWHM'] = '[0, 3, 3]'

        # Add number of bootstraps
        inputs['nBoot'] = 5000

        # Add tau expression
        inputs['tau'] = '1/np.sqrt(nSub)'

        # Create mu1 specification
        mu1 = {}

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

        # Create mu2 specification
        mu2 = {}

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
        fg_nSubs = np.array([100,300,500])

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

        #--------------------------------------------------------------------------------------
        # Save the baseline configuration (Note: This must be the last file output as the bash
        # script used for running simulations on the cluster will take this files existence as
        # a sign to run the next stage of the simulations)
        #--------------------------------------------------------------------------------------
        # Delete fields which vary across simulation
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

        # Save the yml
        with open(os.path.join(simDir,'cfgs','baseline_cfg.yml'), 'w') as outfile:
            yaml.dump(inputs, outfile, default_flow_style=False)

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

        # These are our sample sizes:
        nSubs = np.linspace(40,500,24)

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # New empty inputs directory
        inputs={}

        # Add output directory
        inputs['OutDir'] = OutDir

        # Add number of realizations
        inputs['nReals'] = 1000

        # Add threshold c
        inputs['c'] = 2

        # Add mode
        inputs['mode'] = 2

        # Add range of p-values
        inputs['p'] = 'np.linspace(0,1,21)'

        # Add simulation number (helps identify simulation)
        inputs['simNo'] = simNo

        # Add FWHM
        inputs['FWHM'] = '[0, 3, 3]'

        # Add number of bootstraps
        inputs['nBoot'] = 10000

        # Add tau expression
        inputs['tau'] = '1/np.sqrt(nSub)'

        # Create mu1 specification
        mu1 = {}

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

        # Create mu2 specification
        mu2 = {}

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
        fg_nSubs = np.array([100,300,500])

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

        #--------------------------------------------------------------------------------------
        # Save the baseline configuration (Note: This must be the last file output as the bash
        # script used for running simulations on the cluster will take this files existence as
        # a sign to run the next stage of the simulations)
        #--------------------------------------------------------------------------------------
        # Delete fields which vary across simulation
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

        # Save the yml
        with open(os.path.join(simDir,'cfgs','baseline_cfg.yml'), 'w') as outfile:
            yaml.dump(inputs, outfile, default_flow_style=False)

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

        # These are our sample sizes:
        nSubs = np.linspace(40,500,24)

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # New empty inputs directory
        inputs={}

        # Add output directory
        inputs['OutDir'] = OutDir

        # Add number of realizations
        inputs['nReals'] = 500

        # Add threshold c
        inputs['c'] = 2/3

        # Add mode
        inputs['mode'] = 2

        # Add range of p-values
        inputs['p'] = 'np.linspace(0,1,21)'

        # Add simulation number (helps identify simulation)
        inputs['simNo'] = simNo

        # Add FWHM
        inputs['FWHM'] = '[0, 3, 3]'

        # Add number of bootstraps
        inputs['nBoot'] = 5000

        # Add tau expression
        inputs['tau'] = '1/np.sqrt(nSub)'

        # Create mu1 specification
        mu1 = {}

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

        # Create mu2 specification
        mu2 = {}

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
        fg_nSubs = np.array([100,300,500])

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

        #--------------------------------------------------------------------------------------
        # Save the baseline configuration (Note: This must be the last file output as the bash
        # script used for running simulations on the cluster will take this files existence as
        # a sign to run the next stage of the simulations)
        #--------------------------------------------------------------------------------------
        # Delete fields which vary across simulation
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

        # Save the yml
        with open(os.path.join(simDir,'cfgs','baseline_cfg.yml'), 'w') as outfile:
            yaml.dump(inputs, outfile, default_flow_style=False)


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

        # These are our sample sizes:
        nSubs = np.linspace(40,500,24)

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # New empty inputs directory
        inputs={}

        # Add output directory
        inputs['OutDir'] = OutDir

        # Add number of realizations
        inputs['nReals'] = 1000

        # Add threshold c
        inputs['c'] = 2

        # Add mode
        inputs['mode'] = 3

        # Add range of p-values
        inputs['p'] = 'np.linspace(0,1,21)'

        # Add simulation number (helps identify simulation)
        inputs['simNo'] = simNo

        # Add FWHM
        inputs['FWHM'] = '[0, 3, 3]'

        # Add number of bootstraps
        inputs['nBoot'] = 10000

        # Add tau expression
        inputs['tau'] = '1/np.sqrt(nSub)'

        # Create mu1 specification
        mu1 = {}

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

        # Create mu2 specification
        mu2 = {}

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
        fg_nSubs = np.array([100,300,500])

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

        #--------------------------------------------------------------------------------------
        # Save the baseline configuration (Note: This must be the last file output as the bash
        # script used for running simulations on the cluster will take this files existence as
        # a sign to run the next stage of the simulations)
        #--------------------------------------------------------------------------------------
        # Delete fields which vary across simulation
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

        # Save the yml
        with open(os.path.join(simDir,'cfgs','baseline_cfg.yml'), 'w') as outfile:
            yaml.dump(inputs, outfile, default_flow_style=False)

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

        # These are our sample sizes:
        nSubs = np.linspace(40,500,24)

        # These are our circles centers
        centers = np.arange(-20,32,2)

        # New empty inputs directory
        inputs={}

        # Add output directory
        inputs['OutDir'] = OutDir

        # Add number of realizations
        inputs['nReals'] = 500

        # Add threshold c
        inputs['c'] = 2/3

        # Add mode
        inputs['mode'] = 3

        # Add range of p-values
        inputs['p'] = 'np.linspace(0,1,21)'

        # Add simulation number (helps identify simulation)
        inputs['simNo'] = simNo

        # Add FWHM
        inputs['FWHM'] = '[0, 3, 3]'

        # Add number of bootstraps
        inputs['nBoot'] = 5000

        # Add tau expression
        inputs['tau'] = '1/np.sqrt(nSub)'

        # Create mu1 specification
        mu1 = {}

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

        # Create mu2 specification
        mu2 = {}

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
        fg_nSubs = np.array([100,300,500])

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

        #--------------------------------------------------------------------------------------
        # Save the baseline configuration (Note: This must be the last file output as the bash
        # script used for running simulations on the cluster will take this files existence as
        # a sign to run the next stage of the simulations)
        #--------------------------------------------------------------------------------------
        # Delete fields which vary across simulation
        del inputs['mu2']['center'], inputs['figGen'], inputs['cfgId'], inputs['nSub']

        # Save the yml
        with open(os.path.join(simDir,'cfgs','baseline_cfg.yml'), 'w') as outfile:
            yaml.dump(inputs, outfile, default_flow_style=False)

    print(os.path.join(simDir,'cfgs'))
