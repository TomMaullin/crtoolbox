import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# Table of observed coverage based on true and estimated boundaries
obs_true = pd.read_csv('/home/tommaullin/Documents/ConfRes/tmp/trueBdry_intrp.csv',header=None)
obs_est = pd.read_csv('/home/tommaullin/Documents/ConfRes/tmp/estBdry_intrp.csv',header=None)

# Column headers
colhdr = ['cfgID', 'n', 'distance']+['p='+('%.2f' % p) for p in np.linspace(0,1,21)]

# Assign column headers
obs_true.columns=colhdr
obs_est.columns=colhdr

# List of n and p values
n_values = np.unique(obs_est['n'].values)
d_values = np.unique(obs_est['distance'].values)
p_values = np.linspace(0,1,21)

# Loop through all values of n
for n in n_values:

    # Loop through all values of p
    for p in p_values:

        obs_est_n = obs_est[obs_est['n']==n].sort_values('distance')
        obs_true_n = obs_true[obs_true['n']==n].sort_values('distance')

        # Distances
        distances_est_n = obs_est_n[['distance']].values
        p_est_n = obs_est_n[['p='+('%.2f' % p)]].values
        distances_true_n = obs_true_n[['distance']].values
        p_true_n = obs_true_n[['p='+('%.2f' % p)]].values

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
        plt.savefig(os.path.join('/home/tommaullin/Documents/ConfRes/tmp/', 'distance_vs_obsp_truep'+str(np.int(100*p))+'_n'+str(np.int(n))+'.png'))

        # Clear figure
        plt.clf()

# Loop through all values of n
for n in n_values:

    # Loop through all values of p
    for p in p_values:

        obs_est_d = obs_est[obs_est['distance']==d].sort_values('n')
        obs_true_d = obs_true[obs_true['distance']==d].sort_values('n')

        # Distances
        n_est_d = obs_est_d[['n']].values
        p_est_d = obs_est_d[['p='+('%.2f' % p)]].values
        n_true_d = obs_true_d[['n']].values
        p_true_d = obs_true_d[['p='+('%.2f' % p)]].values

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
        plt.savefig(os.path.join('/home/tommaullin/Documents/ConfRes/tmp/', 'n_vs_obsp_truep'+str(np.int(100*p))+'_n'+str(np.int(n))+'.png'))

        # Clear figure
        plt.clf()