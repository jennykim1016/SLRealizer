# ======================================================================                                    
from __future__ import print_function
from scipy.stats import multivariate_normal
#import plotting
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, MaxNLocator
from numpy import linspace
import matplotlib
import scipy
import scipy.optimize as opt
import math
import photutils
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture
from photutils import DAOStarFinder
import constant
import desc.slrealizer
import scipy.stats
from scipy.stats import moment
from scipy import ndimage
import skimage
import skimage.measure
from skimage.measure import moments

# ======================================================================                           

#global variable that controls the size of the plot
x_min = -5.0
x_max = 5.0
y_min = -5.0
y_max = 5.0
distance = 0.01

number_of_rows = int((x_max - x_min)/distance)
number_of_columns = int((y_max - y_min)/distance)
#=======================================================================

"""
The methods in this file will be called by `deblend.py`. 
This file contains methods that will emulate the null-deblender.
Right now, the working deblender is not being supported -- we only support null deblender.
"""

def null_deblend(image2):
    """
    Given a 2d input array, this method generates an output as null deblender would do.
    Returns the array's zeroth moment(flux), first moments, and second moments
    """
    #skimage.measure.moments return four by four matrix, whose entrees are x_column, y_row's moments. i.e) matrix[0][1] is the first moment of x (which is a little counter intuitive)
    moment_matrix = skimage.measure.moments(image2)
    zeroth_moment = moment_matrix[0][0]
    first_moment_x = x_min + (moment_matrix[0][1] / zeroth_moment) * distance
    first_moment_y = y_min + (moment_matrix[1][0] / zeroth_moment) * distance
    # another comment
    moment_matrix = skimage.measure.moments_central(image2, moment_matrix[1][0]/zeroth_moment, moment_matrix[0][1]/zeroth_moment)
    covariance_matrix = [[moment_matrix[0][2], moment_matrix[1][1]], [moment_matrix[1][1], moment_matrix[2][0]]]
    covariance_matrix /= (zeroth_moment)
    # the volume under multivariate is 10,000. Thus, we have to divide with 10000 to get a right value
    covariance_matrix /= 10000
    flux = zeroth_moment / 10000
    return flux, first_moment_x, first_moment_y, covariance_matrix

def null_deblend_plot(flux, first_moment_x, first_moment_y, covariance_matrix):
    """
    Given the zeroth moment, first moment, and the second moment, this method draws the realistic color plot that represents the shape of null-deblended object.
    """
    x, y = np.mgrid[x_min:x_max:distance, y_min:y_max:distance]
    pos = np.dstack((x, y))
    rv = scipy.stats.multivariate_normal([first_moment_x,first_moment_y], covariance_matrix, allow_singular=True) #FIX BUG     
    image = [[0]*number_of_rows for _ in range(number_of_columns)]
    image = image + rv.pdf(pos)*flux
    return image
    
