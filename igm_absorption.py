"""
Code to process and calculate IGM absorption using different models
"""

import numpy as np
from scipy.special import factorial, gamma as func_gamma, gammaincc as func_gammaincc


################################################################################
# constants

#transitions to n=41
linewv = np.array([1215.67, 1025.72, 972.537, 949.743, 937.803, 930.748, 926.226,
                    923.150, 920.963, 919.352, 918.129, 917.181, 916.429, 915.824, 915.329, 914.919, 914.576, 914.286, 914.039, 913.826, 913.641, 913.480, 913.339, 913.215, 913.104, 913.006, 912.918, 912.839, 912.768, 912.703, 912.645, 912.592, 912.543, 912.499, 912.458, 912.420, 912.385, 912.353, 912.324])

#tau_n/tau_alpha for up to n=41.
Rn = np.array([1.00000000e+00, 2.77633136e-01, 1.32485207e-01,
                    7.80473373e-02, 5.15207101e-02, 3.65562130e-02, 
                    2.72721893e-02, 2.11183432e-02, 1.68224852e-02, 
                    1.37159763e-02, 1.13786982e-02, 9.59763314e-03, 
                    8.19526627e-03, 7.07692308e-03, 6.17159763e-03, 
                    5.42840237e-03, 4.80946746e-03, 4.29053254e-03, 
                    3.84911243e-03, 3.47218935e-03, 3.14733728e-03, 
                    2.86568047e-03, 2.61952663e-03, 2.40414201e-03, 
                    2.21183432e-03, 2.04378698e-03, 1.89289941e-03, 
                    1.75798817e-03, 1.63668639e-03, 1.52781065e-03, 
                    1.42899408e-03, 1.33905325e-03, 1.25798817e-03, 
                    1.18343195e-03, 1.11538462e-03, 1.05266272e-03, 
                    9.95266272e-04, 9.42603550e-04, 8.93491124e-04])[:,None]




#################################################################################

def tau_eff(z):
    z_1= 1.2593846589041013
    a_1= 1.405759341921824
    A1= 0.01549433546364817
    z_2= 5.174708845486191
    a_2= 3.6336748056970474
    a_3= 7.587297773671628
    delta= 0.030496534844542546
    A2=2.0448040499201565

    return np.piecewise(z,[z<z_1, z>=z_1], [lambda x:A1*((1.+x))**a_1, lambda z:A2*((1.+z)/(1.+z_2))**(a_2)*(0.5*(1.+((1.+z)/(1.+z_2))**(1./delta)))**(delta*(a_3-a_2))])

def tau_eff_highz(z):
    '''optical depth valid for z>1.26, use if you know your object is high redshift OR that your wavelength range doesnt go below lambda<3206.34'''
    z_2= 5.174708845486191
    a_2= 3.6336748056970474
    a_3= 7.587297773671628
    delta= 0.030496534844542546
    A2=2.0448040499201565
    return A2*((1.+z)/(1.+z_2))**(a_2)*(0.5*(1.+((1.+z)/(1.+z_2))**(1./delta)))**(delta*(a_3-a_2))


## lyman continuum absorption

def tau_lc_igm(zs, l_obs):
    """Calculate the absorption <912 due to optically thin Lyman alpha forest clouds, following Meiksin 2006."""
    lratio = l_obs/911.8
    tau_igm = 0.805*lratio**3. * (1/lratio - 1/(1+zs))
    return tau_igm


    
    
def tau_lc_lls(zs, l_obs, n0=0.15, beta=1.28, gamma=1.94, nterms = 10):
    """Calculate the absorption <912 due to optically thick Lyman Limit Systems"""
    #compute constant terms
    lratio = l_obs/911.8 #ratio of observed wavelength to lyman limit, also z_LLS
    zs_p1 = 1+zs
    e_neg1 = 0.36787944117144233  #1/e
    
    #check if n0, beta, and gamma are the defaults and if so, use precomputed terms. otherwise compute
    if beta==1.28 and gamma==1.94:
        beta_m1 = 0.28
        gamma_p1 = 2.94
        gam_fn = 0.31340445237536
    else:
        gamma_p1 = gamma + 1
        beta_m1 = beta - 1
        gam_fn = func_gamma(2-beta,1)*func_gammaincc(2-beta,1)
        
    sum1 = _first_lls_sum(beta_m1,nterms)
    sum2 = _second_lls_sum(zs_p1, lratio, gamma_p1=gamma_p1, beta_m1=beta_m1, nterms=nterms)

    tau_lls = n0/(4+gamma-3*beta) * (gam_fn - e_neg1 - sum1) * ( zs_p1**(-3*beta_m1+gamma_p1) * lratio**(3*beta_m1) - np.power(lratio,gamma_p1)) - n0*sum2

    #not valid when tau<1
    #tau_lls[tau_lls<0]=0.
    
    return tau_lls

def _first_lls_sum(beta_m1=0.28, nterms=10):
    #if default beta, nterms return precomputed
    if beta_m1 == 0.28 and nterms==10:
        return -1.3219480209375662
    else:
        # sums from n=0
        n = np.arange(0,nterms)
        return np.sum( (beta_m1) / (n - beta_m1) * (-1)**n / factorial(n) )

def _second_lls_sum(zs_p1, lratio, gamma_p1=2.94, beta_m1=0.28, nterms=10):
    # sums from n=1
    n = np.arange(1,1+nterms)[:,None]
    
    #precomputed values for speed if gamma, beta, nterms are the defaults
    if beta_m1 ==0.28 and gamma_p1 == 2.94 and nterms==10:
        gamma_3n_p1 = np.array([ 0.06,  3.06,  6.06,  9.06, 12.06, 15.06, 18.06, 21.06, 24.06,27.06])[:,None]
        term1 = np.array([-6.48148148e+00,  2.65997872e-02, -2.83116547e-03,
         3.46159020e-04, -4.09909024e-05,  4.51444675e-06,
        -4.57762778e-07,  4.27131736e-08, -3.67775581e-09,
         2.93360030e-10])[:,None]
    #otherwise compute!
    else:
        gamma_3n_p1 = 3*n - gamma_p1
        term1 = beta_m1 * (-1)**n / ( (gamma_3n_p1) * (n-beta_m1) * factorial(n) )
    
    #term 2 must be calculated regardless
    term2 = zs_p1**(-gamma_3n_p1) * np.power(lratio, 3*n) - np.power(lratio,gamma_p1)
    
    return np.sum(term1*term2,axis=0)



def calc_transmission(z_s,l_obs,model='full', nlines = 30, z_lls=None):
    """ This function calculates the transmission function for a source at redshift 'z' for an observed wavelength array of l_obs.  

    'z_s': The redshift of the source object

    'l_obs': Wavelength array to calculate transmission over

    'model': default='full'.  Controls whether the full model (valid for all z) is used or the high-redshift (z>1.26) faster model is used.  both agree in the redshift range z>1.26.  The options are 'full', 'highz', or 'optimise'.  'optimise' automatically chooses the best function based on wavelength range.

    'nlines': optional, int, default=30 (up to transition n=32). The number of lines in the lyman series to use.  Default is nlines=30, maximum is 39 (transition n=41).  We do not recommend having fewer than nlines<20, which agrees to within 1% with nlines=39.

    'z_lls': optional, float.  If this is specified, it sets transmission to zero at z_lls instead of using the average absorption from contribution for lls.     

    """
    
    #require that z_lls<=z
    if z_lls is not None:
        if z_lls>z_s:
            raise ValueError('z_lls must be less than or equal to z_s or None')
        
    #check nlines
    if nlines>39:
        raise ValueError('Maximum value for nlines is 39.')
    

    #specify lya model
    if model=='full':
        tau_lya = tau_eff
    elif model=='highz':
        tau_lya = tau_eff_highz
    elif model=='optimise':
        # check the wavelength range to see if any wavelengths are below zbreak, and if not use more efficient highz model
        zbreak = 1.26
        if np.any(l_obs<zbreak):
            tau_lya = tau_eff
        else:
            tau_lya = tau_eff_highz
    


    #lya forest, lyman series
    zlook = np.outer(1./linewv[0:nlines],l_obs)-1
    tau_laf_i = np.zeros_like(zlook)
    # only compute for l_obs<1216*(1+zs)
    mask_lya = zlook<z_s
    tau_laf_i[mask_lya] = tau_lya(zlook[mask_lya])
    tau_laf_i *= Rn[0:nlines]
    tau_lys = np.sum(tau_laf_i,axis=0)

    #IGM contribution to lyman continuum absorption
    mask_lc = l_obs<911.8*(1.+z_s)
    tau_lc = np.zeros_like(l_obs)
    tau_lc[mask_lc] = tau_lc_igm(z_s,l_obs[mask_lc])

    #LLS absorption, calculate average unles z_lls is specified
    tau_lls = np.zeros_like(l_obs)
    if z_lls is None:
        tau_lls[mask_lc] += tau_lc_lls(z_s,l_obs[mask_lc])
    else:
        lls_mask = l_obs<(911.8*(1+z_lls))
        tau_lls[lls_mask] = np.inf
        
        

    #calculate transmission
    
    return np.exp(-(tau_lys+tau_lc+tau_lls))
