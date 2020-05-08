import unittest
import numpy as np

from solver.Face import Face
from solver.TetrahedronCell import TetrahedronCell

class TestFace(unittest.TestCase):

    f1 = Face(np.array([2., 0., 0.]),
              np.array([0., 2., 0.]),
              np.array([0., 0., 2.]))

    f2 = Face(np.array([0., 0., 0.]),
              np.array([2., 0., 0.]),
              np.array([2., 3., 0.]),
              np.array([0., 3., 0.]))

    f3 = Face(np.array([0., -1., 0.]), np.array([0., 1., 0.]), np.array([1., 0., 0.]))
    f4 = Face(np.array([0., -1., 0.]), np.array([0., 0., 1.]), np.array([1., 0., 0.]))
    f5 = Face(np.array([0., 1., 0.]), np.array([0., 0., 1.]), np.array([1., 0., 0.]))

    f6 = Face(np.array([0., -1., 0.]), np.array([0., 0., 1.]), np.array([0., 1., 0.]))

    f7 = Face(np.array([0., -1., 0.]), np.array([0., 1., 0.]), np.array([-1., 0., 0.]))
    f8 = Face(np.array([0., -1., 0.]), np.array([0., 0., 1.]), np.array([-1., 0., 0.]))
    f9 = Face(np.array([0., 1., 0.]), np.array([0., 0., 1.]), np.array([-1., 0., 0.]))

    faces_1 = [f3, f4, f5, f6]
    faces_2 = [f6, f7, f8, f9]
    c1 = TetrahedronCell(*faces_1)
    c2 = TetrahedronCell(*faces_2)

    def test_n_(self):
        self.assertTrue(np.allclose(self.f1.n_, 1.0 / np.sqrt(3.0) * np.array([1., 1., 1.])))

    def test_normal_1(self):
        self.assertTrue(np.allclose(self.f1.normal(np.array([0., 0., 0.])), 1.0 / np.sqrt(3.0) * np.array([1., 1., 1.])))

    def test_normal_2(self):
        self.assertTrue(np.allclose(self.f1.normal(np.array([4., 4., 4.])), -1.0 / np.sqrt(3.0) * np.array([1., 1., 1.])))

    def test_area(self):
        self.assertEqual(self.f2.area, 6.)

    def test_r_(self):
        self.assertTrue(np.allclose(self.f2.r_, np.array([1., 1.5, 0.])))

    def test_set_cell_1(self):
        self.assertEqual(self.f6.right_cell, self.c2)

    def test_set_cell_2(self):
        self.assertEqual(self.f6.left_cell, self.c1)

    def test_is_same_vertices_1(self):
        test_vertices = [np.array([0., 0., 2.]), np.array([0., 2., 0.]), np.array([2., 0., 0.])]
        self.assertTrue(self.f1.is_same_vertices(*test_vertices))

    def test_is_same_vertices_2(self):
        test_vertices = [np.array([0., 0., 2.]), np.array([0., 4., 0.]), np.array([2., 0., 0.])]
        self.assertFalse(self.f1.is_same_vertices(*test_vertices))

    def test_is_same_vertices_3(self):
        test_vertices = [np.array([0., 0., 2.]), np.array([2., 0., 0.])]
        self.assertFalse(self.f1.is_same_vertices(*test_vertices))

if __name__ == '__main__':
    unittest.main()
