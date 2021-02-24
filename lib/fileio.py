import os
import time
import pandas as pd
import numpy as np

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# This file contains all miscellaneous file i/o functions used by the conf 
# sets code. These functions exist so that basic file handling does not take
# too much space in the bulk of the main code.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Tom Maullin (Last edited 10/11/2020)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ============================================================================
#
# This function takes in the data `data` and appends it to the csv file
# `fname`. If the file does not exist already it creates it. A file lock 
# system is also implemented to ensure the file is not edited by multiple
# jobs at the same time.
#
# ----------------------------------------------------------------------------
#
# This function takes the following inputs:
#
# ----------------------------------------------------------------------------
#
# - `fname`: The filename of the file we wish to append data to.
# - `data`: The data we wish to append to the file.
#
# ============================================================================
def append_to_file(fname, data):

    # Check if file is in use
    fileLocked = True
    while fileLocked:

        try:

            # Create lock file, so other jobs know we are writing to this file
            os.open(fname + ".lock", os.O_CREAT|os.O_EXCL|os.O_RDWR)
            fileLocked = False

        except FileExistsError:

        	# File is still locked
            fileLocked = True

    # Check whether the file exists already
    if not os.path.isfile(fname):

        # Output data to file
        pd.DataFrame(data).to_csv(fname, header=False, index=False)
        
    else:

    	# Read in data from file
    	pd.DataFrame(data).to_csv(fname, mode='a', header=False, index=False)

    # Delete lock file, so other jobs know they can now write to the
    # file
    os.remove(fname + ".lock")

# ============================================================================
#
# The below function takes in a string representing a vector and returns the
# vector as an array.
#
# ============================================================================
def str2vec(c):

    c = str(c)
    c = c.replace("'", "")
    c = c.replace('][', '], [').replace('],[', '], [').replace('] [', '], [')
    c = c.replace('[ [', '[[').replace('] ]', ']]')
    cs = c.split(' ')
    cf = ''
    for i in range(0,len(cs)):
        cs[i]=cs[i].replace(',', '')
        cf=cf + cs[i]
        if i < (len(cs)-1):
            cf = cf + ', '
        
    return(eval(cf))