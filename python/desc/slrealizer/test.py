# ======================================================================
import os, unittest
#import null_deblend
import numpy as np
import scipy
import skimage

# ======================================================================

class TestCase(unittest.TestCase):
    
    """
    Tests SLRealizer class.
    Notes
    -----
    Execute these tests with:
        nosetests
    from anywhere in the module, provided you have run
        pip install nose
    """
    
    # ------------------------------------------------------------------

    def setUp(self):
        # setup 2d array with given dimensions (declared as global variables)
        global x_min, x_max, y_min, y_max, distance
        x_min, x_max, y_min, y_max = 0, 10, 0, 10
        distance = 0.01
        x, y = np.mgrid[x_min:x_max:distance, y_min:y_max:distance]
        pos = np.dstack((x, y))
        number_of_rows = int((x_max-x_min)/distance)
        number_of_columns = int((y_max-y_min)/distance)
        image = [[0]*number_of_rows for _ in range(number_of_columns)]
        # return initialized 2d array
        return image

    def generate_gaussian(self):
        image = self.setUp()
        # declare the center of moment and covariance matrix as a global variable, so that test_multivariate can access them
        global x_com, y_com
        global covariance_matrix
        x_com, y_com = 4.0, 6.0
        covariance_matrix = [[0.8, 0.4], [0.4, 0.6]]
        center_of_mass = [x_com, y_com]
        rv = scipy.stats.multivariate_normal([x_com, y_com], covariance_matrix)
        x, y = np.mgrid[x_min:x_max:distance, y_min:y_max:distance]
        pos = np.dstack((x, y))
        global flux
        flux = 2
        # regardless of what flux is, we know that all x_com, y_com, and covariance_matrix should not change
        image = image + rv.pdf(pos)*flux
        return image

    def test_multivariate(self):
        # this method tests whether we can get the same moment values from the array produced in generate_gaussian method
        image = self.generate_gaussian()
        moment_matrix = skimage.measure.moments(image)
        zeroth_moment = moment_matrix[0][0]
        print('Zeroth moment: ', zeroth_moment)
        first_moment_x = moment_matrix[0][1] / zeroth_moment
        first_moment_y = moment_matrix[1][0] / zeroth_moment
        moment_matrix = skimage.measure.moments_central(image, first_moment_x, first_moment_y)
        covariance_matric = [[moment_matrix[2][0], moment_matrix[1][1]], [moment_matrix[1][1], moment_matrix[0][2]]]
        covariance_matric /= zeroth_moment
        covariance_matric /= 10000
        print('Covariance matrix is: ', covariance_matric)
        self.assertAlmostEqual(distance*first_moment_x, x_com, places=3)
        self.assertAlmostEqual(distance*first_moment_y, y_com, places=3)
        self.assertAlmostEqual(covariance_matric[1][1], covariance_matrix[0][0], places=3)
        self.assertAlmostEqual(covariance_matric[1][0], covariance_matrix[1][0], places=3)
        self.assertAlmostEqual(covariance_matric[0][0], covariance_matrix[1][1], places=3)
        self.assertAlmostEqual(covariance_matric[1][0], covariance_matrix[0][1], places=3)
        return

# ======================================================================

if __name__ == '__main__':

    unittest.main()

# ======================================================================