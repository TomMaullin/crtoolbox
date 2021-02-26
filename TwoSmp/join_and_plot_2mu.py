import os
import numpy as np
import yaml
import glob
import pandas as pd
from matplotlib import pyplot as plt
from lib.fileio import *

def joinAndPlot(OutDir, simNo):

    # Get simulation directory
    simDir = os.path.join(OutDir, 'sim'+str(simNo))

    # Make directory to store configuration files
    cfgFiles = glob.glob(os.path.join(simDir,'cfgs','cfg*.yml'))

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
                                                covp_est.reshape(1,n_p)),\
                                                axis=1)

                # Line for table of true boundary results
                tableLine_true = np.concatenate((np.array([[cfgId,nSub,distance]]),\
                                                 covp_true.reshape(1,n_p)),\
                                                 axis=1)

                # Line for table of estimated boundary interpolation assessed results
                tableLine_est_intrp = np.concatenate((np.array([[cfgId,nSub,distance]]),\
                                                      covp_est_intrp.reshape(1,n_p)),\
                                                      axis=1)

                # Line for table of true boundary interpolation assessed results
                tableLine_true_intrp = np.concatenate((np.array([[cfgId,nSub,distance]]),\
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
                # tableLine_time = pd.read_csv(os.path.join(resDir,'computationTime.csv'), header=None, index_col=None)

                # # If this is the first cfg we've looked at, intialize the results tables
                # if first:

                #     # Initialize time table
                #     table_time = pd.DataFrame(tableLine_time)

                # else:

                #     # Append to existing time table
                #     table_time = table_time.append(pd.DataFrame(tableLine_time))

                # ------------------------------------------------------------------
                # Delete files
                # ------------------------------------------------------------------
                # Delete folder for this simulation
                #shutil.rmtree(os.path.join(OutDir, 'sim'+str(simNo), 'cfg' + str(cfgId)))

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
        colhdr = ['cfgID', 'n', 'distance']+['p='+('%.2f' % p) for p in np.linspace(0,1,21)]

        # Assign column headers
        table_true.columns=colhdr
        table_est.columns=colhdr

        # List of n and p values
        n_values = np.unique(table_est['n'].values)
        d_values = np.unique(table_est['distance'].values)
        p_values = np.linspace(0,1,21)

        # Loop through all values of n
        for n in n_values:

            # Loop through all values of p
            for p in p_values:

                table_est_n = table_est[table_est['n']==n].sort_values('distance')
                table_true_n = table_true[table_true['n']==n].sort_values('distance')

                # Distances
                distances_est_n = table_est_n[['distance']].values
                p_est_n = table_est_n[['p='+('%.2f' % p)]].values
                distances_true_n = table_true_n[['distance']].values
                p_true_n = table_true_n[['p='+('%.2f' % p)]].values

                plt.plot(distances_est_n,p_est_n,color="red",label="Estimated boundary")
                plt.plot(distances_true_n,p_true_n,color="blue",label="True boundary")
                plt.hlines(p, np.min(distances_est_n), np.max(distances_est_n),linestyles='dashed',label="Expected")

                # Title
                plt.title("Coverage (" + str(np.int(100*p)) + "% probability, " + str(int(n)) + " subjects)")

                # Axes
                plt.xlabel("Distance between circles")
                plt.ylabel("Observed coverage")

                # Legend
                plt.legend()

                # Save plots
                plt.savefig(os.path.join(fResDir, 'coverage_vs_distance_p'+str(np.int(100*p))+'_n'+str(np.int(n))+'.png'))

                # Clear figure
                plt.clf()

        # Loop through all values of n
        for n in n_values:

            # Loop through all values of p
            for p in p_values:

                table_est_d = table_est[table_est['distance']==d].sort_values('n')
                table_true_d = table_true[table_true['distance']==d].sort_values('n')

                # Distances
                n_est_d = table_est_d[['n']].values
                p_est_d = table_est_d[['p='+('%.2f' % p)]].values
                n_true_d = table_true_d[['n']].values
                p_true_d = table_true_d[['p='+('%.2f' % p)]].values

                plt.plot(n_est_d,p_est_d,color="red",label="Estimated boundary")
                plt.plot(n_true_d,p_true_d,color="blue",label="True boundary")
                plt.hlines(p, np.min(n_est_d), np.max(n_est_d),linestyles='dashed',label="Expected")

                # Title
                plt.title("Coverage (" + str(np.int(100*p)) + "% probability, distance " + str(int(d)) + ")")

                # Axes
                plt.xlabel("Number of subjects")
                plt.ylabel("Observed coverage")

                # Legend
                plt.legend()

                # Save plots
                plt.savefig(os.path.join(fResDir, 'n_vs_obsp_truep'+str(np.int(100*p))+'_n'+str(np.int(n))+'.png'))

                # Clear figure
                plt.clf()