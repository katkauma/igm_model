This repository contains python code for the IGM absorption model discussed in Kauma et al. (in prep). 

igm_absorption.py contains functions to generate transmission as a function of redshift and observed wavelength.

The transmission function can be called with the following model and arguments:
calc_transmission(z_s, l_obs, nlines=30, z_lls=None, lls_kwargs{})
  * z_s: float. is the source object redshift
  * l_obs: array. is the observed wavelength array
  * n_lines: int, optional. is the number of transitions to inlcude in the lyman alpha forest.  Maximum number is 41, I don't recommend going below 20,
  * z_lls: float, optional.  If specified, this turns off the average LLS absorption and instead puts flux to zero at the location z_LLS.

calc_tranmsission returns a transmission array of len(l_obs).


example.ipynb has example plots and usage.


