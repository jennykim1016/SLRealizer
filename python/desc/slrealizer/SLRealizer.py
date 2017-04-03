import desc.slrealizer
import pandas as pd
import matplotlib.pyplot as plt
import pylab
import matplotlib
import math
from __future__ import print_function

class SLRealizer(object):

    scale_factor = 2

    # Maybe make a separate method for it
    # Best OOP practice?

    def __init__(self, catalog=None, observation="../../data/twinkles_observation_history.csv"):
        """
            Reads in a lens sample catalog and observation data.
            We assume lenses are OM10 lenses and observation file is .csv file
        """
        self.catalog = catalog
        self.observation = pd.read_csv(observation,index_col=0).as_matrix()

    # def plot_lens_random_date(self, lensID = 7176527, convolve=False):
    #     plotting.plot_lens_random_date(self, lensID, convolve)

    def plot_lens_random_date(self, lensID=None, convolve=False):
        if lensID is None:
            print('No lens system selected for plotting.')
            return
        import random
        # Keep randomly selecting epochs until we get one that is not
        # in the 'y' filter:
        filter = 'y'
        while filter == 'y':
            randomIndex = random.randint(0, 200)
            filter = self.observation[randomIndex][1]
        # Now visualize the lens system at the epoch defined by the
        # randomIndex:
        desc.slrealizer.draw_model(self.observation[randomIndex],
                                   self.catalog.get_lens(lensID),
                                   convolve)
        return
