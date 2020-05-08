import unittest
import numpy as np

from solver.TetrahedronCell import TetrahedronCell
from solver.Face import Face

class TestTetrahedronCell(unittest.TestCase):

    f1 = Face(np.array([0., 0., 0.]), np.array([1., 0., 0.]), np.array([0., 1., 0.]))
    f2 = Face(np.array([0., 0., 0.]), np.array([1., 0., 0.]), np.array([0., 0., 1.]))
    f3 = Face(np.array([0., 0., 0.]), np.array([0., 1., 0.]), np.array([0., 0., 1.]))
    f4 = Face(np.array([1., 0., 0.]), np.array([0., 1., 0.]), np.array([0., 0., 1.]))
    faces = [f1, f2, f3, f4]
    c1 = TetrahedronCell(*faces)

    def test_r_(self):
        self.assertTrue(np.allclose(self.c1.r_, np.array([0.25, 0.25, 0.25])))

    def test_vol(self):
        self.assertEqual(self.c1.vol, 1.0 / 6.0)

if __name__ == '__main__':
    unittest.main()
