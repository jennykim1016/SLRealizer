#=========================================

import desc.slrealizer
import deblend
import scipy
import scipy.stats

#=========================================

"""
This file can be used to calculate the statistical distances between two images (to be exact, two images arrays)
"""


def chi_square_distance(image1, image2):
    """
    Given two images, calculate the chi square and its p value between the two.
    """

    chi2, p, dof, ex = scipy.stats.chi2_contingency(image1, image2)
    print('chi squared: ', chi2, 'p: ', p)
    return chi2, p


def KL_distance(image1, image2):
    """
    Given two images, calculate the KL divergence between the two
    2d array is not supported, so we have to flatten the array and compare each pixel in the image1 to the corresponding pixel in the image2.
    """

    return scipy.stats.entropy(image1.ravel(), image2.ravel())
