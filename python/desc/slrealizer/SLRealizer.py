#=====================================================

import desc.slrealizer
import pandas as pd
import matplotlib.pyplot as plt
import pylab
import matplotlib
import math
import skimage
import random
import om10

#=====================================================

class SLRealizer(object):

    """
    Contains the constructor and the key methods for SLRealizer module.
    Constructor reads in an OM10 catalog and observation history.
    Generates the toy catalog, plots the lensed system, and deblends sources using OM10 catalog and observation history.
    """

    def __init__(self, catalog=None, observation="../../data/twinkles_observation_history.csv"):
        """
        Reads in a lens sample catalog and observation data.
        We assume lenses are OM10 lenses and observation file is .csv file
        """
        self.catalog = catalog
        self.observation = pd.read_csv(observation,index_col=0).as_matrix()

    def plot_lens_random_date(self, lensID=None, debug=False, convolve=False):
        """
        Given a specific lens, this code plots a lens after choosing a random observation epoch.
        """

        if lensID is None:
            print 'No lens system selected for plotting.'
            return
        # Keep randomly selecting epochs until we get one that is not in the 'y' filter:
        filter = 'y'
        while filter == 'y':
            randomIndex = random.randint(0, 200)
            filter = self.observation[randomIndex][1]
        # Now visualize the lens system at the epoch defined by the randomIndex:
        desc.slrealizer.draw_model(self.observation[randomIndex],
                                   self.catalog.get_lens(lensID),
                                   convolve, debug)
        return

    def make_catalog(self, num_system = 3, save = True):
        """
        Selects the lensed system just as the real LSST will do, and generates a toy catalog

        Parameters
        ----------
        num_system: int
        Number of systems that the user will request

        save: bool
        If true, the catalog will be saved in the data folder.
        """

        print('From the OM10 catalog, I am selecting LSST lenses')
        self.catalog.select_random(maglim=23.3,area=20000.0,IQ=0.75)
        df = pd.DataFrame(columns=['MJD', 'filter', 'RA', 'RA_err', 'DEC', 'DEC_err', 'x', 'x_com_err', 'y', 'y_com_err', 'flux', 'flux_err', 'qxx', 'qxx_err', 'qyy', 'qyy_err', 'qxy', 'qxy_err', 'psf_sigma', 'sky', 'lensid'])
        for i in xrange(num_system):
            randomIndex = random.randint(0, len(self.catalog.sample))
            lensID = self.catalog.sample[randomIndex]['LENSID']
            filter = 'y'
            #  Keep randomly selecting epochs until we get one that is not in the 'y' filter:
            while filter == 'y':
                randomIndex = random.randint(0, 200)
                filter = self.observation[randomIndex][1]
            data = desc.slrealizer.generate_data(self.catalog.get_lens(lensID), self.observation[randomIndex])
            df.loc[len(df)]= data
        if save:
            print('saving the table with the name catalog.csv. Check your data folder (../../../data/)')
            df.to_csv('../../../data/catalog.csv', index=False)

    def deblend(self, lensID=None, null_deblend=True):
        """
        Given a lens system, this method deblends the source and plots the process of deblending.

        Parameters
        ---------
        lensID : int
        OM10 lens ID which can be used to identify the lensed system
        
        null_deblend : bool
        If true, assumes null deblender. Working deblender is currently not being supported
        """

        if lensID is None:
            print('No lens system selected for calculating the statistics')
            return
        if null_deblend is False:
            print('Sorry, working deblender is not being supported.')
            return
        # Keep randomly selecting epochs until we get one that is not in the 'y' filter:
        filter = 'y'
        while filter == 'y':
            randomIndex = random.randint(0, 200)
            filter = self.observation[randomIndex][1]
        image2 = desc.slrealizer.plot_all_objects(self.observation[randomIndex], self.catalog.get_lens(lensID))
        print('##################### PLOTTING ALL SOURCES ##################################')
        desc.slrealizer.show_color_map(image2)
        flux, first_moment_x, first_moment_y, covariance_matrix = desc.slrealizer.null_deblend(image2)
        print('first:', flux, first_moment_x, first_moment_y, covariance_matrix)
        image = desc.slrealizer.null_deblend_plot(flux, first_moment_x, first_moment_y, covariance_matrix)
        print('##################### AFTER NULL DEBLENDING ##################################')
        plt.clf()
        desc.slrealizer.show_color_map(image)
        print('##################### DEBLENDING DEBLENDED OBJECT #################################')
        flux, first_moment_x, first_moment_y, covariance_matrix = desc.slrealizer.null_deblend(image)
        print('first:',flux, first_moment_x, first_moment_y, covariance_matrix)
        image = desc.slrealizer.null_deblend_plot(flux, first_moment_x, first_moment_y, covariance_matrix)
        plt.clf()
        desc.slrealizer.show_color_map(image)
        print('##################### DEBLENDING DEBLENDED OBJECT #################################')
        flux, first_moment_x, first_moment_y, covariance_matrix = desc.slrealizer.null_deblend(image)
        print('first:',flux, first_moment_x, first_moment_y, covariance_matrix)
        image = desc.slrealizer.null_deblend_plot(flux, first_moment_x, first_moment_y, covariance_matrix)
        plt.clf()
        desc.slrealizer.show_color_map(image)
        print('##################### DIFFERENCE BETWEEN TWO IMAGES ##################################')
        desc.slrealizer.show_color_map(image2-image)
