import numpy as np 
import os
import time
import matplotlib.pyplot as plt
from lib.generateData import *
from lib.boundary import *


# Obtain max and min fields (corresponding to union and intersection
# for excursion sets) ("combinations" of mu)
def get_2mu_cmbtns(mu):

    print('Marker')
    print(mu.shape)
    
    # Get max and min fields
    mu_cup = np.maximum(mu[:1,:,:],mu[1:,:,:])
    mu_cap = np.minimum(mu[:1,:,:],mu[1:,:,:])
    
    # print(muMax.shape)
    # print(muMin.shape)

    # # Concateate
    # mu_cmbtns = np.concatenate((muMax,muMin),axis=0)

    # Return result
    return(mu_cup, mu_cap)

# Get the boundaries of:
# - mu1
# - mu2
# - mu1 intersect mu2
# - mu1 union mu2
# - dmu1 intersect dmu2
def get_2mu_bdry_locs(mu,c,image=False):

    # ----------------------------------------------------------------------
    # Mu 1 and Mu 2
    # ----------------------------------------------------------------------

    # Split into mu1 and mu2 (it is easier to handle the coordinates this
    # way and the time saved via broadcasting is neglible for small axes)
    mu1 = mu[:1,:,:]
    mu2 = mu[1:,:,:]

    # ----
    # Mu 1
    # ----

    # Get the excursion set map of mu1
    mu1_bdry_maps = get_bdry_maps(mu1,c)

    print('usvbso')
    print(mu1.shape)
    print(mu1_bdry_maps)

    # Get a map of the inner and outer boundaries across 2 dimensions.
    # This will be used to obtain the intersection and unions of the
    # boundaries
    mu1_bdry_map = (mu1_bdry_maps[1]['top']['inner']+
                    mu1_bdry_maps[1]['bottom']['inner']+
                    mu1_bdry_maps[2]['top']['inner']+
                    mu1_bdry_maps[2]['bottom']['inner']+
                    mu1_bdry_maps[1]['top']['outer']+
                    mu1_bdry_maps[1]['bottom']['outer']+
                    mu1_bdry_maps[2]['top']['outer']+
                    mu1_bdry_maps[2]['bottom']['outer']) > 0

    # Make an image (for testing purposes only)
    if image:

        # Make a plot of mu 1 boundayr
        plt.figure(0)
        plt.imshow(1*mu1_bdry_map[0,:,:])

    # Obtain mu boundary locations
    mu1_bdry_locs = get_bdry_locs(mu1_bdry_maps)

    # ----
    # Mu 2
    # ----

    # Get the excursion set maps of mu1 and mu2
    mu2_bdry_maps = get_bdry_maps(mu2,c)

    # Get a map of the inner and outer boundaries across 2 dimensions.
    # This will be used to obtain the intersection and unions of the
    # boundaries
    mu2_bdry_map = (mu2_bdry_maps[1]['top']['inner']+
                    mu2_bdry_maps[1]['bottom']['inner']+
                    mu2_bdry_maps[2]['top']['inner']+
                    mu2_bdry_maps[2]['bottom']['inner']+
                    mu2_bdry_maps[1]['top']['outer']+
                    mu2_bdry_maps[1]['bottom']['outer']+
                    mu2_bdry_maps[2]['top']['outer']+
                    mu2_bdry_maps[2]['bottom']['outer']) > 0

    # Make an image (for testing purposes only)
    if image:

        # Make a plot of mu 1 boundayr
        plt.figure(1)
        plt.imshow(1*mu2_bdry_map[0,:,:])

    # Obtain mu boundary locations
    mu2_bdry_locs = get_bdry_locs(mu2_bdry_maps)

    # ----------------------------------------------------------------------
    # Intersection and union
    # ----------------------------------------------------------------------

    # Get union and intersection of mu1 and mu2 excursion sets
    mu_cup, mu_cap = get_2mu_cmbtns(mu)

    print(mu.shape)
    print(mu_cap.shape)
    print(mu_cup.shape)

    # -----
    # Union
    # -----

    # Get the excursion set maps of the combinations
    mu_cup_bdry_maps = get_bdry_maps(mu_cup,c)

    # Make an image (for testing purposes only)
    if image:

        # Get a map of the inenr and outer boundaries across 2 dimensions
        # for display purposes
        mu_cup_bdry_map = (mu_cup_bdry_maps[1]['top']['inner']+
                           mu_cup_bdry_maps[1]['bottom']['inner']+
                           mu_cup_bdry_maps[2]['top']['inner']+
                           mu_cup_bdry_maps[2]['bottom']['inner']+
                           mu_cup_bdry_maps[1]['top']['outer']+
                           mu_cup_bdry_maps[1]['bottom']['outer']+
                           mu_cup_bdry_maps[2]['top']['outer']+
                           mu_cup_bdry_maps[2]['bottom']['outer']) > 0

        # Make a plot of mu 1 boundary
        plt.figure(2)
        plt.imshow(1*mu_cup_bdry_map[0,:,:])

    # Obtain mu boundary locations
    mu_cup_bdry_locs = get_bdry_locs(mu_cup_bdry_maps)

    # ------------
    # Intersection
    # ------------

    # Get the excursion set maps of the combinations
    mu_cap_bdry_maps = get_bdry_maps(mu_cap,c)

    # Make an image (for testing purposes only)
    if image:

        # Get a map of the inenr and outer boundaries across 2 dimensions
        # for display purposes
        mu_cap_bdry_map = (mu_cap_bdry_maps[1]['top']['inner']+
                           mu_cap_bdry_maps[1]['bottom']['inner']+
                           mu_cap_bdry_maps[2]['top']['inner']+
                           mu_cap_bdry_maps[2]['bottom']['inner']+
                           mu_cap_bdry_maps[1]['top']['outer']+
                           mu_cap_bdry_maps[1]['bottom']['outer']+
                           mu_cap_bdry_maps[2]['top']['outer']+
                           mu_cap_bdry_maps[2]['bottom']['outer']) > 0

        # Make a plot of mu 1 boundary
        plt.figure(3)
        plt.imshow(1*mu_cap_bdry_map[0,:,:])

    # Obtain mu boundary locations
    mu_cap_bdry_locs = get_bdry_locs(mu_cap_bdry_maps)

    # ----------------------------------------------------------------------
    # Intersections and unions of boundarys
    # ----------------------------------------------------------------------

    # Get the intersection of the two boundaries
    mu_bdry_cap = mu1_bdry_map & mu2_bdry_map

    # Get the union of the two boundaries
    mu_bdry_cup = mu1_bdry_map | mu2_bdry_map 

    # Make an image (for testing purposes only)
    if image:

        # Make a plot of intersection of boundarys
        plt.figure(4)
        plt.imshow(1*mu_bdry_cap[0,:,:])

        # Make a plot of union of boundarys
        plt.figure(5)
        plt.imshow(1*mu_bdry_cup[0,:,:])

    # Obtain mu boundary intersection and union locations
    mu_bdry_cap_locs = np.where(mu_bdry_cap)
    mu_bdry_cup_locs = np.where(mu_bdry_cup)

    # ----------------------------------------------------------------------
    # Clean up and return
    # ----------------------------------------------------------------------
    # Show plot if we made one
    if image:
        plt.show()

    # Return boundary locations
    return(mu1_bdry_locs, mu2_bdry_locs, mu_cap_bdry_locs, mu_cup_bdry_locs, mu_bdry_cap_locs, mu_bdry_cup_locs)


def testfn():

    mu = get_mu({'type': 'circle2D_2mu', 'center': np.array([[-20,0],[20,0]]), 'fwhm': np.array([3,3]), 'r': np.array([25,25]), 'mag': np.array([3,3])}, np.array([1,100,100]))

    mu_cmbtns = get_2mu_cmbtns(mu)

    c = 2

    get_2mu_bdry_locs(mu,mu_cmbtns,c,image=True)