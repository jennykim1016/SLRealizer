"""
Example unit tests for SLRealizer package
"""
import unittest
import om10
import desc.slrealizer

class SLRealizerTestCase(unittest.TestCase):
    def setUp(self):
        self.db = om10.DB()
        self.db.paint(synthetic=True)

    def tearDown(self):
        pass

    def test_init(self):
        realizer = desc.slrealizer.SLRealizer(catalog=self.db,
            observation="data/twinkles_observation_history.csv")
        self.assertGreater(len(realizer.catalog.sample), 0)
        self.assertGreater(len(realizer.observation), 0)

    # def test_failure(self):
    #     self.assertRaises(TypeError, desc.slrealizer.SLRealizer)
    #     foo = desc.slrealizer.SLRealizer(self.message)
    #     self.assertRaises(RuntimeError, foo.run, True)

if __name__ == '__main__':
    unittest.main()
