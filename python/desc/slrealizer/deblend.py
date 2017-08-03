# ======================================================================
from __future__ import print_function
from scipy.stats import multivariate_normal
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

#global variable that controls the size of the plot#
x_min = 0.0
x_max = 10.0
y_min = 0.0
y_max = 10.0
distance = 0.01

number_of_rows = int((x_max - x_min)/distance)
number_of_columns = int((y_max - y_min)/distance)
# =========================================================================

"""                                                                                                                                                                             
Given a specific date and the OM10 catalog, this deblends the sources that are on the catalog.                                                                                   
Assumes null deblender where all the sources are assumed to be observed as single objects.                                                                                        All the sources are assumed to have Gaussian PSFs.                                                                                                                               
"""

def deblend(currObs, currLens, null_deblend=True):
    """
    If the user wants to see the plot drawn by plotting.py in the debug mode, this code draws it.
    Otherwise, it acts like a wrapper method -- this just calls blend_all_objects.
    """
    if null_deblend:
        image = null_deblending(currObs, currLens, debug)
        show_color_map(image)
    else: # For now, SLRealizer wrapper method does not allow the user to pass null_deblend = False 
        pass
        #image = plot_all_objects(currObs, currLens, debug)
        #blend_all_objects(currObs, currLens, debug, image)

def plot_all_objects(currObs, currLens):
    """
    Given a current observation epoch details and a lensed system, this method blends all the light sources, assuming the gaussian PSF. 
    The blending is done by generating a 2d-gaussian and adding the values to the 2d array.
    Then, that 2d array is passed to the deblending method to find the discrete sources.
    """

    filterLens = currObs[1] + '_SDSS_lens'
    lens_mag = currLens[filterLens]
    galaxy_x, galaxy_y, PSF_sigma = 0, 0, currObs[2]
    x, y = np.mgrid[x_min:x_max:distance, y_min:y_max:distance]
    pos = np.dstack((x, y))
    rv = scipy.stats.multivariate_normal([galaxy_x,galaxy_y], [[PSF_sigma*PSF_sigma, 0], [0, PSF_sigma*PSF_sigma]], allow_singular=True)
    image = [[0]*number_of_rows for _ in range(number_of_columns)]
    image = image + rv.pdf(pos)*math.pow(2.5, desc.slrealizer.return_zeropoint() - currLens[filterLens])
    # iterate for the lens
    for i in xrange(currLens['NIMG']):
        mag_ratio = math.pow(2.5, desc.slrealizer.return_zeropoint() -currLens['MAG'][0][i])
        rv = scipy.stats.multivariate_normal([currLens['XIMG'][0][i],currLens['YIMG'][0][i]], [[PSF_sigma*PSF_sigma, 0], [0, PSF_sigma*PSF_sigma]], allow_singular=True)
        image = image + rv.pdf(pos)*mag_ratio #scale
    return image


"""
We can use the code when we are implementing working deblender, but right now, it is not being used.
def blend_all_objects(currObs, currLens, debug, input_image):
    if debug:
        show_color_map(input_image)
    # detect sources using IRAF routine
    segm = photutils.detect_sources(input_image, 0.1, 8)
    if debug:
        plt.imshow(segm.data, origin='lower', interpolation='nearest')
    daofind = DAOStarFinder(fwhm=currObs[2], threshold=0.1)
    sources = daofind(input_image)
    if debug:
        show_source_position(sources, input_image)

def show_source_position(sources, input_image):
    #Only called in debug mode. Draws the 2d image and mark where the sources are detected.
    positions = (sources['xcentroid'], sources['ycentroid'])
    apertures = CircularAperture(positions, r=4.)
    norm = ImageNormalize(stretch=SqrtStretch())
    plt.imshow(input_image, cmap='Greys', origin='lower', norm=norm)
    apertures.plot(color='red', lw=3.0, alpha=1.0)

"""

def show_color_map(input_image):
    """
    Given a 2d array, this draws a color plot that shows the intensity of each pixel.
    """
    cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap', ['black', 'green', 'yellow'], 256)                                                               
    img2 = plt.imshow(input_image, interpolation='nearest', cmap = cmap2, origin='lower', extent=[x_min, x_max, y_min, y_max], aspect = "auto")                                 
    plt.colorbar(img2,cmap=cmap2)                                                                                                                                               
    plt.show()
