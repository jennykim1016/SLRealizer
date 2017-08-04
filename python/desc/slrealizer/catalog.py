#==============================================

import numpy as np
import pandas
import desc.slrealizer
from fractions import Fraction

#==============================================

"""
This file contains methods that help generate the toy catalog for the SLRealizer.
"""


def generate_data(curr_lens, curr_obs):    
    """
    Given the current observation detail and current lensed system, generate a mock catalog.
    """

    ###
    first_moment_err = 0.01 # unit : percentage
    first_moment_err_std = 0.005
    flux_err = 0.01 # unit : percentage                                                                            
    flux_err_std = 0.005
    second_moment_err = 0.05 # unit : percentage
    second_moment_err_std = 0.01
    pos_err = 0.0 # unit : degree
    pos_err_std = Fraction(1, 3) # Assuming that the three sigma is one degree, one sigma is 0.3333 degree 
    ###

    processed_image = desc.slrealizer.plot_all_objects(curr_obs, curr_lens)
    flux, first_moment_x, first_moment_y, covariance_matrix = desc.slrealizer.null_deblend(processed_image)
    flux_err_calc = flux * noissify_data(flux_err, flux_err_std)
    first_moment_x_err_calc = noissify_data(first_moment_err, first_moment_err_std) * first_moment_x
    first_moment_y_err_calc = noissify_data(first_moment_err, first_moment_err_std) * first_moment_y
    RA, RA_err, DEC, DEC_err = return_coordinate(first_moment_x, first_moment_y)
    MJD, filter, PSF_HWHM, sky = curr_obs[0], curr_obs[1], curr_obs[2], curr_obs[3]
    I_xx, I_xy, I_yy = covariance_matrix[0][0], covariance_matrix[0][1], covariance_matrix[1][1]
    I_xx_err_calc = I_xx * noissify_data(second_moment_err, second_moment_err_std)
    I_xy_err_calc = I_xy * noissify_data(second_moment_err, second_moment_err_std)
    I_yy_err_calc = I_yy * noissify_data(second_moment_err, second_moment_err_std)
    lensID = curr_lens[0]['LENSID']
    return np.array([MJD, filter, RA, RA_err, DEC, DEC_err, first_moment_x, first_moment_x_err_calc, first_moment_y, first_moment_y_err_calc, flux, flux_err_calc, I_xx, I_xx_err_calc, I_yy, I_yy_err_calc, I_xy, I_xy_err_calc, PSF_HWHM, sky, lensID])

def noissify_data(mean, stdev):
    """
    Given a mean and a standard deviation of a measurement, add a noise to the data
    """
    return np.random.normal(loc=mean, scale=stdev)

def return_coordinate(first_moment_x, first_moment_y):
    """
    This method returns the RA and DEC values that take the x and y offset into account.
    RA and DEC values are provided in the observational history,
    and this method also assumes a random error in the measurement (+1/-1 deg)
    """
    real_coordinate = desc.slrealizer.return_obs_RA_DEC()
    RA = real_coordinate.ra.deg
    DEC = real_coordinate.dec.deg
    # copied from OM10
    # make sure first_moment_x is com
    RA += first_moment_x/np.cos(DEC)
    DEC += first_moment_y/3600
    #===========================
    pos_err = 0.0 # unit : degree                                                                                                             
    pos_err_std = Fraction(1, 3)
    #===========================
    RA_err = noissify_data(pos_err, pos_err_std)
    DEC_err = noissify_data(pos_err, pos_err_std)
    return RA, RA_err, DEC, DEC_err
