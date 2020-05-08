import unittest
import numpy as np

from solver.Face import Face

class TestFace(unittest.TestCase):

    def test_normal_vector(self):
        f = Face(np.array([0., 0., 0.]),
                 np.array([2., 0., 0.]),
                 np.array([0., 3., 0.]))
        self.assertTrue(np.allclose(f.n_, np.array([0., 0., 1.])))

    def test_surface_area(self):
        f = Face(np.array([0., 0., 0.]),
                 np.array([2., 0., 0.]),
                 np.array([2., 3., 0.]),
                 np.array([0., 3., 0.]))
        self.assertEqual(f.area, 6.)

    def test_centroid(self):
        f = Face(np.array([0., 0., 0.]),
                 np.array([2., 0., 0.]),
                 np.array([2., 3., 0.]),
                 np.array([0., 3., 0.]))
        self.assertTrue(np.allclose(f.c_, np.array([1., 1.5, 0.])))

if __name__ == '__main__':
    unittest.main()
