import numpy as np 

# Powerset function, taken from:
# https://stackoverflow.com/questions/1482308/how-to-get-all-subsets-of-a-set-powerset
def powerset(s):

    # Get length of x
    x = len(s)

    # Masks is a list of `aliases' for the elements of s. I.e. for the
    # i^th element of s, the corresponding element of masks is 2^(i+1).
    masks = [1 << i for i in range(x)]
    
    # a << x is the same as a*(2^x) so we are looping over length of
    # powerset here
    for i in range(1,1 << x):

        # Yeild makes a generator obect
        # zip(masks,s) is pairing each element of s with a power of 2
        yield np.array([ss for mask, ss in zip(masks, s) if i & mask])
        
