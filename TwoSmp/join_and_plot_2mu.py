import os
import numpy as np
import yaml
import glob
import shutil
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
    # Simulations 1-8: Circle/Squares moving closer
    #
    # --------------------------------------------------------------------------
    #
    # In this simulation setting, we are interested in moving two circles (or 
    # squares) of equal diameter close to one another. For this reason, we vary
    # the circles (squares) center and, as usual, the number of subjects.
    #
    # ==========================================================================
    if simNo in [1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16]:

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
                tableLine_time = pd.read_csv(os.path.join(resDir,'times.csv'), header=None, index_col=None)

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
        append_to_file(os.path.join(fResDir,'times.csv'), table_time)

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
        table_true_intrp.columns=colhdr
        table_est_intrp.columns=colhdr

        # List of n and p values
        n_values = np.unique(table_est_intrp['n'].values)
        d_values = np.unique(table_est_intrp['distance'].values)
        p_values = np.linspace(0,1,21)

        # Loop through all values of n
        for n in n_values:

            # Loop through all values of p
            for p in p_values:

                table_est_n = table_est_intrp[table_est_intrp['n']==n].sort_values('distance')
                table_true_n = table_true_intrp[table_true_intrp['n']==n].sort_values('distance')

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
                if simNo in [1,2]:
                    plt.xlabel("Distance between circles")
                else:
                    plt.xlabel("Distance between squares")
                plt.ylabel("Observed coverage")

                # Make axis a bit clearer
                plt.ylim((np.min(p_true_n)-0.02,1))
                
                # Legend
                plt.legend()

                # Save plots
                plt.savefig(os.path.join(fResDir, 'd_vs_obsp_truep'+str(np.int(100*p))+'_n'+str(np.int(n))+'.png'))

                # Clear figure
                plt.clf()

        # Loop through all values of distance
        for d in d_values:

            # Loop through all values of p
            for p in p_values:

                table_est_d = table_est_intrp[table_est_intrp['distance']==d].sort_values('n')
                table_true_d = table_true_intrp[table_true_intrp['distance']==d].sort_values('n')

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
                
                # Make axis a bit clearer
                plt.ylim((np.min(p_true_d)-0.02,1))
                
                # Legend
                plt.legend()

                # Save plots
                plt.savefig(os.path.join(fResDir, 'n_vs_obsp_truep'+str(np.int(100*p))+'_d'+str(np.int(d))+'.png'))

                # Clear figure
                plt.clf()

    if simNo in [9,10]:

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
                # Get number of subjects and fwhm for second noise field
                # ------------------------------------------------------------------

                # Number of subjects
                nSub = inputs['nSub']

                # FWHM for second noise field
                FWHM2 = str2vec(inputs['FWHM2'])[1]

                # ------------------------------------------------------------------
                # Add coverage probabilities to table
                # ------------------------------------------------------------------
                # Line for table of estimated boundary results
                tableLine_est = np.concatenate((np.array([[cfgId,nSub,FWHM2]]),\
                                                covp_est.reshape(1,n_p)),\
                                                axis=1)

                # Line for table of true boundary results
                tableLine_true = np.concatenate((np.array([[cfgId,nSub,FWHM2]]),\
                                                 covp_true.reshape(1,n_p)),\
                                                 axis=1)

                # Line for table of estimated boundary interpolation assessed results
                tableLine_est_intrp = np.concatenate((np.array([[cfgId,nSub,FWHM2]]),\
                                                      covp_est_intrp.reshape(1,n_p)),\
                                                      axis=1)

                # Line for table of true boundary interpolation assessed results
                tableLine_true_intrp = np.concatenate((np.array([[cfgId,nSub,FWHM2]]),\
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
                tableLine_time = pd.read_csv(os.path.join(resDir,'times.csv'), header=None, index_col=None)

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
        append_to_file(os.path.join(fResDir,'times.csv'), table_time)

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
        colhdr = ['cfgID', 'n', 'FWHM2']+['p='+('%.2f' % p) for p in np.linspace(0,1,21)]

        # Assign column headers
        table_true_intrp.columns=colhdr
        table_est_intrp.columns=colhdr

        # List of n and p values
        n_values = np.unique(table_est_intrp['n'].values)
        f_values = np.unique(table_est_intrp['FWHM2'].values)
        p_values = np.linspace(0,1,21)

        # Loop through all values of n
        for n in n_values:

            # Loop through all values of p
            for p in p_values:

                table_est_n = table_est_intrp[table_est_intrp['n']==n].sort_values('FWHM2')
                table_true_n = table_true_intrp[table_true_intrp['n']==n].sort_values('FWHM2')

                # FWHMs
                fwhm2_est_n = table_est_n[['FWHM2']].values
                p_est_n = table_est_n[['p='+('%.2f' % p)]].values
                fwhm2_true_n = table_true_n[['FWHM2']].values
                p_true_n = table_true_n[['p='+('%.2f' % p)]].values

                plt.plot(fwhm2_est_n,p_est_n,color="red",label="Estimated boundary")
                plt.plot(fwhm2_true_n,p_true_n,color="blue",label="True boundary")
                plt.hlines(p, np.min(fwhm2_est_n), np.max(fwhm2_est_n),linestyles='dashed',label="Expected")

                # Title
                plt.title("Coverage (" + str(np.int(100*p)) + "% probability, " + str(int(n)) + " subjects)")

                # Axes
                plt.xlabel("FWHM for 2nd noise field")
                plt.ylabel("Observed coverage")

                # Make axis a bit clearer
                plt.ylim((np.min(p_true_n)-0.02,1))
                
                # Legend
                plt.legend()

                # Save plots
                plt.savefig(os.path.join(fResDir, 'fwhm_vs_obsp_truep'+str(np.int(100*p))+'_n'+str(np.int(n))+'.png'))

                # Clear figure
                plt.clf()

        # Loop through all values of fwhm2
        for f in f_values:

            # Loop through all values of p
            for p in p_values:

                table_est_f = table_est_intrp[table_est_intrp['FWHM2']==f].sort_values('n')
                table_true_f = table_true_intrp[table_true_intrp['FWHM2']==f].sort_values('n')

                # n and p for this fwhm
                n_est_f = table_est_f[['n']].values
                p_est_f = table_est_f[['p='+('%.2f' % p)]].values
                n_true_f = table_true_f[['n']].values
                p_true_f = table_true_f[['p='+('%.2f' % p)]].values

                plt.plot(n_est_f,p_est_f,color="red",label="Estimated boundary")
                plt.plot(n_true_f,p_true_f,color="blue",label="True boundary")
                plt.hlines(p, np.min(n_est_f), np.max(n_est_f),linestyles='dashed',label="Expected")

                # Title
                plt.title("Coverage (" + str(np.int(100*p)) + "% probability, FWHM2 " + ('%.2f' % f) + ")")

                # Axes
                plt.xlabel("Number of subjects")
                plt.ylabel("Observed coverage")
                
                # Make axis a bit clearer
                plt.ylim((np.min(p_true_f)-0.02,1))
                
                # Legend
                plt.legend()

                # Save plots
                plt.savefig(os.path.join(fResDir, 'n_vs_obsp_truep'+str(np.int(100*p))+'_FWHM'+('%.2f' % f)+'.png'))

                # Clear figure
                plt.clf()

    # Simulated over correlation range
    if simNo in [17,18]:

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

                # Correlation between noise fields
                corr = np.float(inputs['noiseCorr'])

                # ------------------------------------------------------------------
                # Add coverage probabilities to table
                # ------------------------------------------------------------------
                # Line for table of estimated boundary results
                tableLine_est = np.concatenate((np.array([[cfgId,nSub,corr]]),\
                                                covp_est.reshape(1,n_p)),\
                                                axis=1)

                # Line for table of true boundary results
                tableLine_true = np.concatenate((np.array([[cfgId,nSub,corr]]),\
                                                 covp_true.reshape(1,n_p)),\
                                                 axis=1)

                # Line for table of estimated boundary interpolation assessed results
                tableLine_est_intrp = np.concatenate((np.array([[cfgId,nSub,corr]]),\
                                                      covp_est_intrp.reshape(1,n_p)),\
                                                      axis=1)

                # Line for table of true boundary interpolation assessed results
                tableLine_true_intrp = np.concatenate((np.array([[cfgId,nSub,corr]]),\
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
                tableLine_time = pd.read_csv(os.path.join(resDir,'times.csv'), header=None, index_col=None)

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
        append_to_file(os.path.join(fResDir,'times.csv'), table_time)

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
        colhdr = ['cfgID', 'n', 'corr']+['p='+('%.2f' % p) for p in np.linspace(0,1,21)]

        # Assign column headers
        table_true_intrp.columns=colhdr
        table_est_intrp.columns=colhdr

        # List of n and p values
        n_values = np.unique(table_est_intrp['n'].values)
        c_values = np.unique(table_est_intrp['corr'].values)
        p_values = np.linspace(0,1,21)

        # Loop through all values of n
        for n in n_values:

            # Loop through all values of p
            for p in p_values:

                table_est_n = table_est_intrp[table_est_intrp['n']==n].sort_values('corr')
                table_true_n = table_true_intrp[table_true_intrp['n']==n].sort_values('corr')

                # Correlations
                corr_est_n = table_est_n[['corr']].values
                p_est_n = table_est_n[['p='+('%.2f' % p)]].values
                corr_true_n = table_true_n[['corr']].values
                p_true_n = table_true_n[['p='+('%.2f' % p)]].values

                print(corr_est_n,p_est_n)

                plt.plot(corr_est_n,p_est_n,color="red",label="Estimated boundary")
                plt.plot(corr_true_n,p_true_n,color="blue",label="True boundary")
                plt.hlines(p, np.min(corr_est_n), np.max(corr_est_n),linestyles='dashed',label="Expected")

                # Title
                plt.title("Coverage (" + str(np.int(100*p)) + "% probability, " + str(int(n)) + " subjects)")

                # Axes
                plt.xlabel("Correlation between noise fields")
                plt.ylabel("Observed coverage")

                # Make axis a bit clearer
                plt.ylim((np.min(p_true_n)-0.02,1))
                
                # Legend
                plt.legend()

                # Save plots
                plt.savefig(os.path.join(fResDir, 'corr_vs_obsp_truep'+str(np.int(100*p))+'_n'+str(np.int(n))+'.png'))

                # Clear figure
                plt.clf()

        # Loop through all values of correlation
        for c in c_values:

            # Loop through all values of p
            for p in p_values:

                table_est_c = table_est_intrp[table_est_intrp['corr']==c].sort_values('n')
                table_true_c = table_true_intrp[table_true_intrp['corr']==c].sort_values('n')

                # n and p for this correlation
                n_est_c = table_est_c[['n']].values
                p_est_c = table_est_c[['p='+('%.2f' % p)]].values
                n_true_c = table_true_c[['n']].values
                p_true_c = table_true_c[['p='+('%.2f' % p)]].values

                plt.plot(n_est_c,p_est_c,color="red",label="Estimated boundary")
                plt.plot(n_true_c,p_true_c,color="blue",label="True boundary")
                plt.hlines(p, np.min(n_est_c), np.max(n_est_c),linestyles='dashed',label="Expected")

                # Title
                plt.title("Coverage (" + str(np.int(100*p)) + "% probability, noise correlation " + ('%.2f' % c) + ")")

                # Axes
                plt.xlabel("Number of subjects")
                plt.ylabel("Observed coverage")
                
                # Make axis a bit clearer
                plt.ylim((np.min(p_true_c)-0.02,1))
                
                # Legend
                plt.legend()

                # Save plots
                plt.savefig(os.path.join(fResDir, 'n_vs_obsp_truep'+str(np.int(100*p))+'_corr'+('%.2f' % c)+'.png'))

                # Clear figure
                plt.clf()



    # Simulated over gradient range
    if simNo in [19,20]:

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
                # Get number of subjects and gradient of ramps
                # ------------------------------------------------------------------

                # Number of subjects
                nSub = inputs['nSub']

                # Gradient of ramps
                grad = (np.float(inputs['mu1']['b']) - np.float(inputs['mu1']['a']))/2

                # ------------------------------------------------------------------
                # Add coverage probabilities to table
                # ------------------------------------------------------------------
                # Line for table of estimated boundary results
                tableLine_est = np.concatenate((np.array([[cfgId,nSub,grad]]),\
                                                covp_est.reshape(1,n_p)),\
                                                axis=1)

                # Line for table of true boundary results
                tableLine_true = np.concatenate((np.array([[cfgId,nSub,grad]]),\
                                                 covp_true.reshape(1,n_p)),\
                                                 axis=1)

                # Line for table of estimated boundary interpolation assessed results
                tableLine_est_intrp = np.concatenate((np.array([[cfgId,nSub,grad]]),\
                                                      covp_est_intrp.reshape(1,n_p)),\
                                                      axis=1)

                # Line for table of true boundary interpolation assessed results
                tableLine_true_intrp = np.concatenate((np.array([[cfgId,nSub,grad]]),\
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
                tableLine_time = pd.read_csv(os.path.join(resDir,'times.csv'), header=None, index_col=None)

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
        append_to_file(os.path.join(fResDir,'times.csv'), table_time)

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
        colhdr = ['cfgID', 'n', 'grad']+['p='+('%.2f' % p) for p in np.linspace(0,1,21)]

        # Assign column headers
        table_true_intrp.columns=colhdr
        table_est_intrp.columns=colhdr

        # List of n and p values
        n_values = np.unique(table_est_intrp['n'].values)
        g_values = np.unique(table_est_intrp['grad'].values)
        p_values = np.linspace(0,1,21)

        # Loop through all values of n
        for n in n_values:

            # Loop through all values of p
            for p in p_values:

                table_est_n = table_est_intrp[table_est_intrp['n']==n].sort_values('grad')
                table_true_n = table_true_intrp[table_true_intrp['n']==n].sort_values('grad')

                # Covariances
                grad_est_n = table_est_n[['grad']].values
                p_est_n = table_est_n[['p='+('%.2f' % p)]].values
                grad_true_n = table_true_n[['grad']].values
                p_true_n = table_true_n[['p='+('%.2f' % p)]].values

                print(grad_est_n,p_est_n)

                plt.plot(grad_est_n,p_est_n,color="red",label="Estimated boundary")
                plt.plot(grad_true_n,p_true_n,color="blue",label="True boundary")
                plt.hlines(p, np.min(grad_est_n), np.max(grad_est_n),linestyles='dashed',label="Expected")

                # Title
                plt.title("Coverage (" + str(np.int(100*p)) + "% probability, " + str(int(n)) + " subjects)")

                # Axes
                plt.xlabel("Slope (per 50 voxels)")
                plt.ylabel("Observed coverage")

                # Make axis a bit clearer
                plt.ylim((np.min(p_true_n)-0.02,1))
                
                # Legend
                plt.legend()

                # Save plots
                plt.savefig(os.path.join(fResDir, 'grad_vs_obsp_truep'+str(np.int(100*p))+'_n'+str(np.int(n))+'.png'))

                # Clear figure
                plt.clf()

        # Loop through all values of fwhm2
        for g in g_values:

            # Loop through all values of p
            for p in p_values:

                table_est_g = table_est_intrp[table_est_intrp['grad']==g].sort_values('n')
                table_true_g = table_true_intrp[table_true_intrp['grad']==g].sort_values('n')

                # n and p for this gradient
                n_est_g = table_est_g[['n']].values
                p_est_g = table_est_g[['p='+('%.2f' % p)]].values
                n_true_g = table_true_g[['n']].values
                p_true_g = table_true_g[['p='+('%.2f' % p)]].values

                plt.plot(n_est_g,p_est_g,color="red",label="Estimated boundary")
                plt.plot(n_true_g,p_true_g,color="blue",label="True boundary")
                plt.hlines(p, np.min(n_est_g), np.max(n_est_g),linestyles='dashed',label="Expected")

                # Title
                plt.title("Coverage (" + str(np.int(100*p)) + "% probability, gradient " + ('%.2f' % g) + ")")

                # Axes
                plt.xlabel("Number of subjects")
                plt.ylabel("Observed coverage")
                
                # Make axis a bit clearer
                plt.ylim((np.min(p_true_g)-0.02,1))
                
                # Legend
                plt.legend()

                # Save plots
                plt.savefig(os.path.join(fResDir, 'n_vs_obsp_truep'+str(np.int(100*p))+'_grad'+('%.2f' % g)+'.png'))

                # Clear figure
                plt.clf()

    # Simulated over noise magnitude range
    if simNo in [21,22]:

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
                # Get number of subjects and magnitude of second noise field
                # ------------------------------------------------------------------

                # Number of subjects
                nSub = inputs['nSub']

                # Magnitude of second noise field
                mag = np.float(inputs['noise2']['mag'])

                # ------------------------------------------------------------------
                # Add coverage probabilities to table
                # ------------------------------------------------------------------
                # Line for table of estimated boundary results
                tableLine_est = np.concatenate((np.array([[cfgId,nSub,mag]]),\
                                                covp_est.reshape(1,n_p)),\
                                                axis=1)

                # Line for table of true boundary results
                tableLine_true = np.concatenate((np.array([[cfgId,nSub,mag]]),\
                                                 covp_true.reshape(1,n_p)),\
                                                 axis=1)

                # Line for table of estimated boundary interpolation assessed results
                tableLine_est_intrp = np.concatenate((np.array([[cfgId,nSub,mag]]),\
                                                      covp_est_intrp.reshape(1,n_p)),\
                                                      axis=1)

                # Line for table of true boundary interpolation assessed results
                tableLine_true_intrp = np.concatenate((np.array([[cfgId,nSub,mag]]),\
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
                tableLine_time = pd.read_csv(os.path.join(resDir,'times.csv'), header=None, index_col=None)

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

        # ----------------------------------------------------------------------
        # Make figures
        # ----------------------------------------------------------------------

        # Column headers
        colhdr = ['cfgID', 'n', 'mag']+['p='+('%.2f' % p) for p in np.linspace(0,1,21)]

        # Assign column headers
        table_true_intrp.columns=colhdr
        table_est_intrp.columns=colhdr

        # List of n and p values
        n_values = np.unique(table_est_intrp['n'].values)
        m_values = np.unique(table_est_intrp['mag'].values)
        p_values = np.linspace(0,1,21)

        # Loop through all values of n
        for n in n_values:

            # Loop through all values of p
            for p in p_values:

                table_est_n = table_est_intrp[table_est_intrp['n']==n].sort_values('mag')
                table_true_n = table_true_intrp[table_true_intrp['n']==n].sort_values('mag')

                # Covariances
                mag_est_n = table_est_n[['mag']].values
                p_est_n = table_est_n[['p='+('%.2f' % p)]].values
                mag_true_n = table_true_n[['mag']].values
                p_true_n = table_true_n[['p='+('%.2f' % p)]].values

                print(mag_est_n,p_est_n)

                plt.plot(mag_est_n,p_est_n,color="red",label="Estimated boundary")
                plt.plot(mag_true_n,p_true_n,color="blue",label="True boundary")
                plt.hlines(p, np.min(mag_est_n), np.max(mag_est_n),linestyles='dashed',label="Expected")

                # Title
                plt.title("Coverage (" + str(np.int(100*p)) + "% probability, " + str(int(n)) + " subjects)")

                # Axes
                plt.xlabel("Magitude of second noise field")
                plt.ylabel("Observed coverage")

                # Make axis a bit clearer
                plt.ylim((np.min(p_true_n)-0.02,1))
                
                # Legend
                plt.legend()

                # Save plots
                plt.savefig(os.path.join(fResDir, 'mag_vs_obsp_truep'+str(np.int(100*p))+'_n'+str(np.int(n))+'.png'))

                # Clear figure
                plt.clf()

        # Loop through all values of noise magnitude
        for m in m_values:

            # Loop through all values of p
            for p in p_values:

                table_est_m = table_est_intrp[table_est_intrp['mag']==m].sort_values('n')
                table_true_m = table_true_intrp[table_true_intrp['mag']==m].sort_values('n')

                # n and p for this gradient
                n_est_m = table_est_m[['n']].values
                p_est_m = table_est_m[['p='+('%.2f' % p)]].values
                n_true_m = table_true_m[['n']].values
                p_true_m = table_true_m[['p='+('%.2f' % p)]].values

                plt.plot(n_est_m,p_est_m,color="red",label="Estimated boundary")
                plt.plot(n_true_m,p_true_m,color="blue",label="True boundary")
                plt.hlines(p, np.min(n_est_m), np.max(n_est_m),linestyles='dashed',label="Expected")

                # Title
                plt.title("Coverage (" + str(np.int(100*p)) + "% probability, noise magnitude " + ('%.2f' % m) + ")")

                # Axes
                plt.xlabel("Number of subjects")
                plt.ylabel("Observed coverage")
                
                # Make axis a bit clearer
                plt.ylim((np.min(p_true_m)-0.02,1))
                
                # Legend
                plt.legend()

                # Save plots
                plt.savefig(os.path.join(fResDir, 'n_vs_obsp_truep'+str(np.int(100*p))+'_mag'+('%.2f' % m)+'.png'))

                # Clear figure
                plt.clf()