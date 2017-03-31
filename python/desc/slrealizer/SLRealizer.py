import desc.slrealizer
import pandas as pd
import matplotlib.pyplot as plt
import pylab
import matplotlib
import math

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
            print 'No lens system selected for plotting.'
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

"""
    def determineAlpha(self, mag_ratio):
        if(mag_ratio>1):
            quasar_alpha = 1/mag_ratio
            lens_alpha = 1
        else:
            quasar_alpha = 1
            lens_alpha = mag_ratio
        return quasar_alpha, lens_alpha

    def determineScale(self, lensX, sourceX, lensY, sourceY):
        if (abs(lensX[0]+sourceX))>(abs(sourceX)):
            plotX = abs(lensX[0]+sourceX)
        else:
            plotX = abs(sourceX)
        if (abs(lensY[0]+sourceY))>(abs(sourceY)):
            plotY = abs(lensY[0]+sourceY)
        else:
            plotY = abs(sourceY)
        return plotX, plotY

    def plotFigureOnMatplotlib(self, currObs, convolve, quasar_alpha, lens_alpha, sourceX, sourceY, lensX, lensY, plotX, plotY):
        fig = plt.figure(figsize=(4, 16))
        scale_factor = 2
        for i in range(4):
            if(currObs[1]=='g'):
                circleColor = 'g'
            else:
                circleColor = 'r'
            plotNumStr = "4"+"1"+str(i+1)
            plotNum = int(plotNumStr)
            fig.subplots_adjust(top=1.5)
            sub = fig.add_subplot(plotNum)
            #sub.set_ylim(-plotY-2*currObs[2]/scale_factor, plotY+2*currObs[2]/scale_factor)
            sub.set_ylim(-3, 3)
            #sub.set_xlim(-plotX-2*currObs[2]/scale_factor, plotX+2*currObs[2]/scale_factor)
            sub.set_xlim(-3, 3)
            sub.set_xlabel('xPosition')
            sub.set_ylabel('yPosition')
            source = plt.Circle((sourceX, sourceY), radius=currObs[2]/scale_factor, alpha=quasar_alpha, fc=circleColor, linewidth=0)
            if(convolve):
                # For now, we just let lensFWHM to be zero.
                #lensFWHM = currLens['REFF']
                lensFWHM = 0.0
                lens = plt.Circle((lensX[i]+sourceX, lensY[i]+sourceY), radius=convolve(currObs[2], lensFWHM)/scale_factor, alpha=lens_alpha, fc=circleColor, linewidth=0)
            else:
                lens = plt.Circle((lensX[i]+sourceX, lensY[i]+sourceY), radius=currObs[2]/scale_factor, alpha=lens_alpha, fc=circleColor, linewidth=0)
            fig.gca().add_patch(lens)
            fig.gca().add_patch(source)
            sub.set_title('Observation ' + str(i+1) + ' with ' + 'filter ' + currObs[1] + ' on MJD ' + str(currObs[0]))
            #seeing = plt.Circle(((-plotY-2*currObs[2])/scale_factor, (plotX-2*currObs[2])/scale_factor), radius=currObs[2]/scale_factor, alpha=0.1, fc='k')
            seeing = plt.Circle((-2.5, -2.5), radius=currObs[2]/scale_factor, alpha=0.1, fc='k')
            sub.legend((source, seeing, lens),('source', 'seeing', 'lens'))
            plt.gca().add_patch(seeing)


    def drawPlot(self, currObs, lensID = 7176527, convolve = False):
        currLens = self.lens.get_lens(lensID)
        #obsHist has MJD Filter FWHM 5sigmag
        filterQuasar = currObs[1] + '_SDSS_quasar'
        filterLens = currObs[1] + '_SDSS_lens'
        lens_mag = currLens[filterLens]
        quasar_mag = currLens[filterQuasar]
        mag_ratio = math.pow(2.5, -lens_mag+quasar_mag)
        quasar_alpha, lens_alpha = self.determineAlpha(mag_ratio)
        sourceX = currLens['XSRC'][0]
        sourceY = currLens['YSRC'][0]
        lensX = currLens['XIMG'][0]
        lensY = currLens['YIMG'][0]
        plotX, plotY = self.determineScale(lensX, sourceX, lensY, sourceY)
        self.plotFigureOnMatplotlib(currObs, convolve, quasar_alpha, lens_alpha, sourceX, sourceY, lensX, lensY, plotX, plotY)

    def convolve(obsFWHM, initialFWHM = 0.0):
        seeing = obsFWHM
        seeingsigma = seeingToSig(seeing)
        initSigma = seeingToSig(initialFWHM)
        convolveSigma = seeingsigma + initSigma
        return convolveSigma

    def seeingToSig(self, seeing):
        return seeing / np.sqrt(8 * np.log(2))

    def sigmaTofwhm(self, sigma):
        return sigma * np.sqrt(8 * np.log(2))
"""
