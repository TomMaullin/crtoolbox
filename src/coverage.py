import os
import sys
import time
import numpy as np
from lib.boundary import *
from lib.fileio import *
from bootstrap import *
from lib.set_theory import powerset

"""
Function to check for inclusion violations of the estimated boundary. This
function takes in FcHat^+, FcHat^-, the original data, means, threshold, 
tau and the quantile levels, a. It returns an array of binary variables
indicating whether or not a violation was observed for each quantile.

Inputs:
    FcHat_plus: Inner confidence region
    FcHat_minus: Outer confidence region
    muhats: Estimated means
    sigmas: Estimated standard deviations
    mus: Original means
    n: Number of subjects
    c: Threshold
    a: Quantile levels
    m: Number of conditions
    X: Design matrix (optional)
    L: Contrast matrix (optional)

Outputs:
    success: Binary variable indicating whether or not a violation was 
             observed for each quantile (using either the binary maps or
             intepolation).
    dict_success: Dictionary of binary variables indicating whether or not
                  a violation was observed for each quantile. The keys are:
     - binary_success: Binary variable indicating whether or not a violation
                       was observed based solely on the binary maps.
     - interp_success: Binary variable indicating whether or not a violation
                       was observed based solely on interpolation.
"""
def check_violations(FcHat_plus, FcHat_minus, muhats, sigmas, mus, n, c, a, m=1, X=None, L=None):

    # -------------------------------------------------------------------
    # Check if a is an array
    # -------------------------------------------------------------------

    # If a is not an array
    if not isinstance(a, np.ndarray):
            
        # Make it one
        a = np.array([a])

    # -------------------------------------------------------------------
    # Load in data if necessary
    # -------------------------------------------------------------------

    # If muhats is a string, assume it is a filename and load it in
    if isinstance(muhats, str):
            
        # Load in muhats
        muhats = read_image(muhats)

    # If sigmas is a string, assume it is a filename and load it in
    if isinstance(sigmas, str):
                
        # Load in sigmas
        sigmas = read_image(sigmas)

    # If mus is a string, assume it is a filename and load it in
    if isinstance(mus, str):

        # Load in mus
        mus = read_image(mus)

    # If FcHat_plus is a string, assume it is a filename and load it in
    if isinstance(FcHat_plus, str):
        
        # Load in FcHat_plus
        FcHat_plus = read_image(FcHat_plus)

    # If FcHat_minus is a string, assume it is a filename and load it in    
    if isinstance(FcHat_minus, str):
            
        # Load in FcHat_minus
        FcHat_minus = read_image(FcHat_minus)

    # -------------------------------------------------------------------
    # Reshape data if necessary
    # -------------------------------------------------------------------

    # If the first dimension of muhats is not m
    if muhats.shape[0] != m:

        # Reshape muhats so that the leading dimension is 1
        muhats = muhats.reshape((1,) + muhats.shape)

    # If the first dimension of sigmas is not m
    if sigmas.shape[0] != m:
            
        # Reshape sigmas so that the leading dimension is 1
        sigmas = sigmas.reshape((1,) + sigmas.shape)

    # If mus has fewer dimensions than muhats
    if mus.ndim < muhats.ndim:

        # Reshape mus so that the leading dimension is 1
        mus = mus.reshape((1,) + mus.shape)

    # If the first dimension of FcHat_plus is not len(a)
    if FcHat_plus.shape[0] != len(a):

        # Reshape FcHat_plus so that the leading dimension is len(a)
        FcHat_plus = FcHat_plus.reshape((len(a),) + FcHat_plus.shape)
    
    # If the first dimension of FcHat_minus is not len(a)
    if FcHat_minus.shape[0] != len(a):
            
        # Reshape FcHat_minus so that the leading dimension is len(a)
        FcHat_minus = FcHat_minus.reshape((len(a),) + FcHat_minus.shape)

    # Get number of dimensions
    D = mus.ndim - 1

    # If X is none, make it a vector of ones
    if X is None:
        X = np.ones((n,1))

    # If L is none, make it a vector of ones
    if L is None:
        L = np.ones((1,1))

    # Work out tau
    tau = np.sqrt(L.T @ np.linalg.pinv(X.T @ X) @ L)[0,0]

    # -------------------------------------------------------------------
    # Get minimum field
    # -------------------------------------------------------------------
    # This is named cap as the excursion set of the minimum field is
    # the intersection of all fields (\cap in latex)
    cap_mu = np.amin(mus,axis=0)

    # Obtain Fc
    Fc = cap_mu > c

    # -------------------------------------------------------------------
    # Some set logic to work out violations
    # -------------------------------------------------------------------

    # Obtain FcHat^+\FcHat based on the estimated boundary. This variable
    # has axes corresponding to [pvalue, field dimensions]
    FcHatp_sub_Fc = FcHat_plus & ~Fc

    # Obtain Fc\FcHat^- based on the estimated boundary. This variable
    # has axes corresponding to [pvalue, field dimensions]
    Fc_sub_FcHatm = Fc & ~FcHat_minus

    # Get a tuple of the form (1,...,D)
    D_tuple = tuple(np.arange(D)+1)

    # Record if we saw a violation in the estimated boundary based sets
    binary_success = 1-(np.any(FcHatp_sub_Fc,axis=D_tuple) | np.any(Fc_sub_FcHatm,axis=D_tuple))

    # -------------------------------------------------------------------
    # Get stat along the Fc boundary
    # -------------------------------------------------------------------

    # Obtain statistic, g
    g = ((muhats-c)/(sigmas*tau))

    # -------------------------------------------------------------------
    # Boundary locations for FcHat
    # -------------------------------------------------------------------

    # Get boolean map for the boundary of Fc
    Fc_bdry_map = get_bdry_maps(cap_mu, c)

    # Get coordinates for the boundary of Fc
    Fc_bdry_locs = get_bdry_locs(Fc_bdry_map)

    # Empty dict to store mu along true boundary
    mu_dFc = {}

    # Loop through fields
    for i in (np.arange(m)+1):

        # -------------------------------------------------------------------
        # Mu along dFc
        # -------------------------------------------------------------------

        # Obtain Mu along Fc
        mu_dFc[str(i)] = get_bdry_values_concat(mus[i-1,...], Fc_bdry_locs)

    # -------------------------------------------------------------------
    # Boundary partitions 
    # -------------------------------------------------------------------

    # Get list of possible alphas to be considered
    alphas=list(powerset(np.arange(m)+1))

    # Locations of the dalpha boundaries
    dalphaFc_locs = {}

    # Loop through alphas
    for alpha in alphas:

        # Loop through possible elements of alpha
        for i in (np.arange(m)+1):

            # Get indices representing where mui and muhati are less
            # than or equal to c if i is in alpha.
            if i in alpha:

                # Mu
                in_dalphaFc_i = (mu_dFc[str(i)][:,1] <= c)

            # Get indices representing where mui is greater than c
            # if i is not in alpha.
            else:

                # Mu
                in_dalphaFc_i = (mu_dFc[str(i)][:,1] > c)

            # Get a boolean index telling us whether each location in dFc
            # belongs to d^\alpha Fc.
            if i == 1:

                # Initial running product for mu
                in_dalphaFc = np.array(in_dalphaFc_i)

            else:

                # Update running product for mu
                in_dalphaFc = in_dalphaFc*in_dalphaFc_i

        # Convert boolean 1/0s into coordinates
        dalphaFc_loc = np.where(in_dalphaFc)[0]

        # Save locations
        dalphaFc_locs[np.array2string(alpha)] = dalphaFc_loc



    # -------------------------------------------------------------------
    # Get residuals, mu and muhat along boundary partitions
    # -------------------------------------------------------------------

    # Empty dict for mu
    mu_dFc_partitioned = {}

    # Loop through boundary partitions
    for alpha in alphas:

        # New empty dict for mu on dalpha Fc
        mu_dalphaFc = {}

        # Get dalpha locations
        dalphaFc_loc = dalphaFc_locs[np.array2string(alpha)]

        # Loop through i in alpha getting values for interpolation
        for i in alpha:

            # ------------------------------------------------------
            # Mu on true boundary; dFc
            # ------------------------------------------------------

            # Get mu for field i along dFc
            mui_dFc = mu_dFc[str(i)]

            # Save residuals for field i along dalphaFc
            mu_dalphaFc[str(i)] = mui_dFc[dalphaFc_loc,:]


        # Save mu and muhat for alpha
        mu_dFc_partitioned[np.array2string(alpha)] = mu_dalphaFc


    # -------------------------------------------------------------------
    # Get weights from mu and muhat along boundary partitions for 
    # interpolation
    # -------------------------------------------------------------------

    # Empty dicts for weights
    weights_dFc = {}

    # Loop through boundary partitions
    for alpha in alphas:

        # New empty dicts for weights
        weights_dalphaFc = {}

        # Loop through i in alpha getting values for interpolation
        for i in alpha:

            # Get mui and muhati along dalpha Fc
            mui_dalphaFc = mu_dFc_partitioned[np.array2string(alpha)][str(i)]

            # For true boundary
            dalphaFc_mui_weights = get_bdry_weights_concat(mui_dalphaFc, c)

            # Save weights
            weights_dalphaFc[str(i)] = dalphaFc_mui_weights

        # Save weights
        weights_dFc[np.array2string(alpha)] = weights_dalphaFc

    # -------------------------------------------------------------------
    # Get stat along the Fc boundary
    # -------------------------------------------------------------------

    g_dFc = {}

    # Loop through fields
    for i in (np.arange(m)+1):

        # Get the values for gi along dFc
        g_dFc[str(i)] = get_bdry_values_concat(g[i-1,...], Fc_bdry_locs)
        
        # -------------------------------------------------------------------
        # Mu along dFc
        # -------------------------------------------------------------------

        # Obtain Mu along Fc
        mu_dFc[str(i)] = get_bdry_values_concat(mus[i-1,...], Fc_bdry_locs)

    # Boolean to tell us if this is the first alpha we've looked at
    Firstalpha = True

    # Interpolate g
    for alpha in alphas:

        # Boolean to tell us if this is the first value of i we're looking at
        Firsti = True

        # Check if we have recorded values for this boundary
        if np.array2string(alpha) in dalphaFc_locs:

            # Loop through all fields relevent to boundary.
            for i in alpha:

                # Get gi on dFc
                gi_dFc = g_dFc[str(i)] 

                # Get locations for dalphaFc
                dalphaFc_loc = dalphaFc_locs[np.array2string(alpha)]

                # Get gi on dalphaFc
                gi_dalphaFc = gi_dFc[dalphaFc_loc,:]

                # Get weights for interpolation
                dalphaFc_weightsi = weights_dFc[np.array2string(alpha)][str(i)]

                # Interpolate gi
                gi_dalphaFc_interp = get_bdry_vals_interpolated_concat(gi_dalphaFc,dalphaFc_weightsi)

                # If this is the first i we've looked at in alpha
                if Firsti:

                    # set the minimum interpolated values to gi_dalphaFc_interp
                    ming_dalphaFc_interp = gi_dalphaFc_interp 

                    # Set first to false
                    Firsti = False

                else:

                    # Take the minimum between current gi_dalphaFc_interp and previous
                    ming_dalphaFc_interp = np.minimum(ming_dalphaFc_interp, gi_dalphaFc_interp) 

            # -------------------------------------------------------------------
            # Get stat along the Fc boundary
            # -------------------------------------------------------------------

            # If this is the first alpha we've looked at, set the stat array to ming_dalphaFc_interp
            if Firstalpha:

                stat_dFc = ming_dalphaFc_interp

                # Now we've seen an alpha
                Firstalpha = False

            # Otherwise concatenate it on
            else:

                stat_dFc = np.concatenate((ming_dalphaFc_interp,stat_dFc),axis=-1)

    # -------------------------------------------------------------------
    # Check whether there were any boundary violations using interpolated
    # boundary values (checking if voxels had values corresponding to no
    # violations, etc)
    # -------------------------------------------------------------------

    # Reshaping for binary comparison
    a = a.reshape(np.prod(a.shape),1)
    stat_dFc2 = stat_dFc.reshape(1,np.prod(stat_dFc.shape))

    # Perform lower check on stat map using thresholds based on the
    # estimated boundary
    bdry_lowerCheck = stat_dFc2 >= -a

    # Perform upper check on stat map using thresholds based on the
    # estimated boundary
    bdry_upperCheck = stat_dFc2 <= a

    # -------------------------------------------------------------------
    # Work out whether we observed successful sets.
    # -------------------------------------------------------------------
    # Record if we saw a violation via the interpolation method
    interp_success = np.all(bdry_lowerCheck,axis=(1)) & np.all(bdry_upperCheck,axis=(1))

    # -------------------------------------------------------------------
    # Success as a whole
    # -------------------------------------------------------------------
    # Final success
    success = interp_success & binary_success

    # Return result
    return(success, {'Interp: ': interp_success, 'Binary: ': binary_success})
