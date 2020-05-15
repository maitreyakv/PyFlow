import unittest
import numpy as np

from solver.Node import Node
from solver.Face import Face
from solver.TetrahedronCell import TetrahedronCell

class TestFace(unittest.TestCase):

    n1 = Node(2., 0., 0., 101)
    n2 = Node(0., 2., 0., 102)
    n3 = Node(0., 0., 2., 103)
    n4 = Node(0., 0., 0., 104)
    n5 = Node(2., 3., 0., 105)
    n6 = Node(0., 3., 0., 106)
    n7 = Node(0., -1., 0., 107)
    n8 = Node(0., 1., 0., 108)
    n9 = Node(1., 0., 0., 109)
    n10 = Node(0., 0., 1., 110)
    n11 = Node(-1., 0., 0., 111)
    n12 = Node(2., 0., 0., 112)

    f1 = Face(n1, n2, n3)

    f2 = Face(n4, n1, n5, n6)

    f3 = Face(n7, n8, n9)
    f4 = Face(n7, n10, n9)
    f5 = Face(n8, n10, n9)

    f6 = Face(n7, n10, n8)

    f7 = Face(n7, n8, n11)
    f8 = Face(n7, n10, n11)
    f9 = Face(n8, n10, n11)

    faces_1 = [f3, f4, f5, f6]
    nodes_1 = [n7, n8, n9, n10]
    faces_2 = [f6, f7, f8, f9]
    nodes_2 = [n7, n8, n10, n11]
    c1 = TetrahedronCell(faces_1, nodes_1)
    c2 = TetrahedronCell(faces_2, nodes_2)

    f10 = Face(n2, n12, n3)

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

    def test_eq_1(self):
        self.assertTrue(self.f1 == self.f10)

    def test_eq_2(self):
        self.assertFalse(self.f1 == self.f2)

    def test_eq_3(self):
        self.assertFalse(self.f3 == self.f4)

    def test_hash_1(self):
        self.assertTrue(hash(self.f1) == hash(self.f10))

    def test_hash_2(self):
        self.assertFalse(hash(self.f1) == hash(self.f2))

    def test_hash_3(self):
        self.assertFalse(hash(self.f3) == hash(self.f4))

    def test_other_cell_1(self):
        self.assertEqual(self.f6.other_cell(self.c1), self.c2)

    def test_other_cell_2(self):
        self.assertEqual(self.f6.other_cell(self.c2), self.c1)

if __name__ == '__main__':
    unittest.main()
