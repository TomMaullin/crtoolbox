import os
import numpy as np
import yaml
import glob
import pandas as pd
from lib.fileio import *

def joinAndPlot(OutDir, simNo):

    # Get simulation directory
    simDir = os.path.join(OutDir, 'sim'+str(simNo))

    # Make directory to store configuration files
    cfgFiles = glob.glob(os.path.join(simDir,'cfgs','*.yml'))

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
            covp_est = np.mean(obs_est,axis=0)[:]

            # Get the coverage probabilities from the observed results for the
            # true boundary
            covp_true = np.mean(obs_true,axis=0)[:]

            # Get the coverage probabilities from the observed results for the
            # estimated boundary (for coverage assessed using interpolation)
            covp_est_intrp = np.mean(obs_est_intrp,axis=0)[:]

            # Get the coverage probabilities from the observed results for the
            # true boundary (for coverage assessed using interpolation)
            covp_true_intrp = np.mean(obs_true_intrp,axis=0)[:]

            # ------------------------------------------------------------------
            # Get number of subjects and distance between radii
            # ------------------------------------------------------------------

            # Number of subjects
            nSub = inputs['nSub']

            # Distance between circle centers
            distance = np.sum(eval(inputs['mu2']['center'])-eval(inputs['mu1']['center']))

            # ------------------------------------------------------------------
            # Add coverage probabilities to table
            # ------------------------------------------------------------------
            # Line for table of estimated boundary results
            tableLine_est = np.concatenate((np.array([[cfgId,nSub,distance]]),\
                                            covp_est.reshape(1,10)),\
                                            axis=1)

            # Line for table of true boundary results
            tableLine_true = np.concatenate((np.array([[cfgId,nSub,distance]]),\
                                             covp_true.reshape(1,10)),\
                                             axis=1)

            # Line for table of estimated boundary interpolation assessed results
            tableLine_est_intrp = np.concatenate((np.array([[cfgId,nSub,distance]]),\
                                                  covp_est_intrp.reshape(1,10)),\
                                                  axis=1)

            # Line for table of true boundary interpolation assessed results
            tableLine_true_intrp = np.concatenate((np.array([[cfgId,nSub,distance]]),\
                                                   covp_true_intrp.reshape(1,10)),\
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

        # ----------------------------------------------------------------------
        # Sort and save to csv
        # ----------------------------------------------------------------------
        # Make final results results directory
        fResDir = os.path.join(OutDir, 'sim'+str(simNo), 'FinalResults')
        if not os.path.exists(fResDir):
            os.mkdir(fResDir)

        # Save times table
        append_to_file(os.path.join(fResDir,'times.csv'), table_time)

        # Save estimated boundary results table
        append_to_file(os.path.join(fResDir,'estBdry.csv'), table_est)

        # Save true boundary results table
        append_to_file(os.path.join(fResDir,'trueBdry.csv'), table_true)

        # Save estimated boundary (with interpolation) results table
        append_to_file(os.path.join(fResDir,'estBdry_intrp.csv'), table_est_intrp)

        # Save true boundary (with interpolation) results table
        append_to_file(os.path.join(fResDir,'trueBdry_intrp.csv'), table_true_intrp)

