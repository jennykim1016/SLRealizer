import pandas as pd
import matplotlib.pyplot as plt
import pylab
import matplotlib
import math
import numpy as np

def draw_model(currObs, currLens, convolved=False):
    """
        Given data of current observation, current lens, this method draws plot of lensed systems.
        If convolve is True, this method convolves two gaussian distributions - lens galaxy's PSF assuming circular model and seeing of the observation date - when plotting the lensed system.
    """
    #obsHist has MJD Filter FWHM 5sigmag
    filterQuasar = currObs[1] + '_SDSS_quasar'
    filterLens = currObs[1] + '_SDSS_lens'
    # Position of lensing galaxy - scalar
    lensX = currLens['XSRC'][0]
    lensY = currLens['YSRC'][0]
    # Position of lensed image - double array(4)
    sourceX = currLens['XIMG'][0]
    sourceY = currLens['YIMG'][0]
    if (currObs[1]=='u'):
        circleColor = 'violet'
    elif (currObs[1]=='g'):
        circleColor = 'blue'
    elif (currObs[1]=='r'):
        circleColor = 'green'
    elif (currObs[1]=='i'):
        circleColor = 'orange'
    elif (currObs[1]=='z'):
        circleColor = 'red'
    else:
        raise ValueError('Unknown filter name '+currObs[1])

    plt.axis('scaled')
    plt.ylim(-3, 3)
    plt.xlim(-3, 3)
    plt.xlabel('xPosition')
    plt.ylabel('yPosition')
    plt.title('Observation with ' + 'filter ' + currObs[1] + ' on MJD ' + str(currObs[0]))
    
    scale_factor = 2
    PSF_HWHM = currObs[2]/scale_factor

    for i in range(4):
        
        print i
        print currLens['MAG'][0]
        print currLens['NIMG']
        #print "In 'draw_model', mag_ratio, quasar_alpha, lens_alpha =", \
        #    mag_ratio, quasar_alpha, lens_alpha
        lens_mag = currLens[filterLens]
        quasar_mag = currLens[filterQuasar] + 2.5*np.log10(currLens['MAG'][0][i])
        print 'Quasar Mag, Lens Mag'
        print quasar_mag, lens_mag
        
        mag_ratio = math.pow(2.5, -lens_mag+quasar_mag)
        # assume that lens_alpha is always greater than quasar_alpha
        # lens_alpha is always one(1)
        quasar_alpha, lens_alpha = determine_alpha(mag_ratio)
        print 'Quasar Alpha, Lens Alpha'
        print quasar_alpha, lens_alpha
        # no convolution by default
        galaxy_HWHM = currObs[2]
        PSF_HWHM = currObs[2]
        if(convolved):
            lens_FWHM = currLens['REFF']
            galaxy_HWHM = convolve(PSF_HWHM, lens_FWHM)/scale_factor
            # Bug: the LENS position is in the CENTER of the field.
            #      The SOURCE position is UNOBSERVABLE.
            #      The IMAGE positions are in teh OM10 catalog in the
            #      XIMG and YIMG columns.
        # Draw lens galaxy:
        print 'Quasar Position X, Y', \
            lensX+sourceX[i], lensY+sourceY[i]
        source = plt.Circle((lensX+sourceX[i], lensY+sourceY[i]),
                            radius=PSF_HWHM,
                            alpha=quasar_alpha,
                            fc=circleColor, linewidth=0)
        plt.gca().add_patch(source)
        #seeing = plt.Circle(((-plotY-2*currObs[2])/scale_factor, (plotX-2*currObs[2])/scale_factor), radius=PSF_HWHM, alpha=0.1, fc='k')
        seeing = plt.Circle((-2.5, -2.5),
                            radius=PSF_HWHM,
                            alpha=0.1,
                            fc='gray')
        plt.gca().add_patch(seeing)
        lens = plt.Circle((lensX, lensY),
                  radius=galaxy_HWHM,
                  alpha=lens_alpha,
                  fc='black', linewidth=3)
        print 'Lens Position X, Y', \
            lensX, lensY
        plt.gca().add_patch(lens)
        plt.legend((source, seeing, lens),
           ('QSO images', 'PSF', 'Lens galaxy'), fontsize=10)


    #plotX, plotY = determine_scale(lensX, sourceX, lensY, sourceY)
    # In order to improve function readabilty, actual plotting involving matplotlib is decomposed
    #plot_figure_on_matplotlib(currObs, convolve, sourceX, sourceY, lensX, lensY, plotX, plotY)



def plot_figure_on_matplotlib(currObs, convolve, sourceX, sourceY, lensX, lensY, plotX, plotY):
    """
        Actual plotting involving matplotlib happens here.
        Based on the filter's wavelength, we chose five different colors to represent the lensed system.
    """
    #fig = plt.figure(1)


def determine_alpha(flux_ratio):
    """
        Given flux ratio, return alpha value for each quasar and lens.
        Brighter object will have alpha = 1, and fainter object will have magnitude = 1/flux_ratio
    """
    """
    Bug: No docstring
    Bug: Each image should have its own alpha, because the different
         images have different magnifications (these are in the OM10 "MAG" columns, and need to be applied to each image's sourceX
         magnitude.) =>? just to make sure, different magnifications but same instrument?
    Bug: function name should be "determine_alpha" to be PEP8 compliant
    """
    if(flux_ratio>1):
        quasar_alpha = 1/flux_ratio
        lens_alpha = 1
    else:
        quasar_alpha = 1
        lens_alpha = flux_ratio
    return quasar_alpha, lens_alpha

"""
def determine_scale(lensX, sourceX, lensY, sourceY):
    if (abs(lensX[0]+sourceX))>(abs(sourceX)):
        plotX = abs(lensX[0]+sourceX)
    else:
        plotX = abs(sourceX)
    if (abs(lensY[0]+sourceY))>(abs(sourceY)):
        plotY = abs(lensY[0]+sourceY)
    else:
        plotY = abs(sourceY)
    return plotX, plotY
"""

def convolve(obsFWHM, initialFWHM=0.0):
    """
        Assume Gaussian Distributions, convolve two gaussians and return the new gaussian distribution's standard deviation
    """
    fwhm_sigma = fwhm_to_sig(obsFWHM)
    initSigma = fwhm_to_sig(initialFWHM)
    convolveSigma = fwhm_sigma + initSigma
    return convolveSigma

def fwhm_to_sig(fwhm):
    """
        Given full width half maximum(fwhm), return sigma(standard deviation)
    """
    return fwhm / np.sqrt(8 * np.log(2))

def sigma_to_fwhm(sigma):
    """
        Given sigma(standard deviation), return full width half maximum(fwhm)
        """
    return sigma * np.sqrt(8 * np.log(2))
