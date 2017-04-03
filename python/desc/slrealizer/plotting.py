import pandas as pd
import matplotlib.pyplot as plt
import pylab
import matplotlib
import math
import numpy as np
np.set_printoptions(threshold='nan')

def draw_model(currObs, currLens, convolved=False):
    """
        Given data of current observation, current lens, this method draws plot of lensed systems.
        If convolve is True, this method convolves two gaussian distributions - lens galaxy's PSF assuming circular model and seeing of the observation date - when plotting the lensed system.
    """
    #obsHist has columns : MJD Filter FWHM 5sigma
    filterQuasar = currObs[1] + '_SDSS_quasar'
    filterLens = currObs[1] + '_SDSS_lens'
    # Position of lensing galaxy - scalar
    lensX = currLens['XSRC'][0]
    lensY = currLens['YSRC'][0]
    # Position of lensed image - double array(4)
    sourceX = currLens['XIMG'][0]
    sourceY = currLens['YIMG'][0]
    
    # Choose color to represent lens system
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

    # some tricks that make plots easy-to-read
    plt.axis('scaled')
    plt.ylim(-3, 3)
    plt.xlim(-3, 3)
    plt.xlabel('xPosition')
    plt.ylabel('yPosition')
    plt.title('Observation with ' + 'filter ' + currObs[1] + ' on MJD ' + str(currObs[0]))

    # a magic number that scales lenses, quasars, and PSF circles to have appropriate size
    scale_factor = 2
    PSF_HWHM = currObs[2]/scale_factor

    # There are two or four quasar images in OM10 catalog
    for i in range(4):
        #print "In 'draw_model', mag_ratio, quasar_alpha, lens_alpha =", \
        #    mag_ratio, quasar_alpha, lens_alpha
        
        lens_mag = currLens[filterLens]
        quasar_mag = currLens[filterQuasar]
        if currLens['MAG'][0][i] is not 0:
            quasar_mag = currLens[filterQuasar] - 2.5*np.log10(abs(currLens['MAG'][0][i]))
        flux_ratio = math.pow(2.5, -lens_mag+quasar_mag)
        # We adjust alpha values based on the flux ratio of images and lens
        quasar_alpha, lens_alpha = determine_alpha(flux_ratio)


        # no convolution by default
        galaxy_HWHM = currObs[2]
        PSF_HWHM = currObs[2]
        if(convolved):
            lens_FWHM = currLens['REFF']
            galaxy_HWHM = convolve(PSF_HWHM, lens_FWHM)/scale_factor

        source = plt.Circle((lensX+sourceX[i], lensY+sourceY[i]),
                            radius=PSF_HWHM,
                            alpha=quasar_alpha,
                            fc=circleColor, linewidth=0)
        plt.gca().add_patch(source)


        if i is 0:
            seeing = plt.Circle((-2.5, -2.5),
                                radius=PSF_HWHM,
                                alpha=0.1,
                                fc='gray')
            plt.gca().add_patch(seeing)
            lens = plt.Circle((lensX, lensY),
                              radius=galaxy_HWHM,
                              alpha=lens_alpha,
                              fc=circleColor, hatch = '/')
            plt.gca().add_patch(lens)
            plt.legend((source, seeing, lens),
                       ('QSO images', 'PSF', 'Lens galaxy'), fontsize=10)


def determine_alpha(flux_ratio):
    """
        Given flux ratio, return alpha value for each quasar and lens.
        Brighter object will have alpha = 1, and fainter object will have magnitude of 1/flux_ratio
    """
    if(flux_ratio>1):
        quasar_alpha = 1/flux_ratio
        lens_alpha = 1
    else:
        quasar_alpha = 1
        lens_alpha = flux_ratio
    return quasar_alpha, lens_alpha

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
