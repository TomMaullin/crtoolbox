
# ----------------------------------------------------------------------------------------------
# Test your understanding: The images we have looked at so far can all be found in the 
# `./demo/data/example1` folder. Another set of CRs can be found in `./demo/data/example2` 
# folder under similar filenames. Using the functions described above, have a look at these 
# files. If these were your analysis results, which CRs do you think would be preferable to 
# observe? And why?
# ----------------------------------------------------------------------------------------------

# Load in pre-computed confidence regions.
Upper_CR_fname = os.path.join(os.getcwd(),'demo','data','example2','Upper_CR_Example2.nii.gz')
Lower_CR_fname = os.path.join(os.getcwd(),'demo','data','example2','Lower_CR_Example2.nii.gz')

# Load in a pre-computed point estimate of Ac.
estimated_Ac_fname = os.path.join(os.getcwd(),'demo','data','example2','Estimated_Ac_Example2.nii.gz')

# Load in a pre-computed mask.
mask_fname = os.path.join(os.getcwd(),'demo','data','mask.nii.gz')

# Display confidence regions in interactive slice plot.
display_crs(estimated_Ac_fname, Upper_CR_fname, Lower_CR_fname, mask_fname, mode='Sagittal')







# ----------------------------------------------------------------------------------------------
# Choose your noise and signal
# ----------------------------------------------------------------------------------------------
my_noise = Noise(var=8)
my_signal = CircleSignal(r=20, mag=3)







# ----------------------------------------------------------------------------------------------
# Let's perform a regression to get the images we need. You can use the same `regression`
# function as in the 2D example to perform this step.
# ----------------------------------------------------------------------------------------------
betahat_files, std_beta_files, resid_files = regression(y_files, X, out_dir)






# ----------------------------------------------------------------------------------------------
# Have a look at these CRs by writing some code in the box below.
# ----------------------------------------------------------------------------------------------
display_crs(estimated_ac_file, upper_cr_file, lower_cr_file, mask_file)





