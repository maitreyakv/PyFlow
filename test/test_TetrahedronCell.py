import unittest
import numpy as np

from solver.Node import Node
from solver.Face import Face
from solver.TetrahedronCell import TetrahedronCell

class TestTetrahedronCell(unittest.TestCase):

    n1 = Node(0., 0., 0., 101)
    n2 = Node(1., 0., 0., 102)
    n3 = Node(0., 1., 0., 103)
    n4 = Node(0., 0., 1., 104)
    n5 = Node(0., 0., 1., 105)

    n6 = Node(999., 0., 0., 106)
    n7 = Node(1000., 0., 0., 107)
    n8 = Node(999., 1., 0., 108)
    n9 = Node(999., 0., 1., 109)

    f1 = Face(n1, n2, n3)
    f2 = Face(n1, n2, n4)
    f3 = Face(n1, n3, n4)
    f4 = Face(n2, n3, n4)
    f5 = Face(n2, n3, n5)
    faces1 = [f1, f2, f3, f4]
    nodes1 = [n1, n2, n3, n4]
    faces2 = [f1, f2, f3, f5]
    nodes2 = [n1, n2, n3, n5]
    c1 = TetrahedronCell(faces1, nodes1)
    c2 = TetrahedronCell(faces2[::-1], nodes2[::-1])

    f6 = Face(n6, n7, n8)
    f7 = Face(n6, n7, n9)
    f8 = Face(n6, n8, n9)
    f9 = Face(n7, n8, n9)
    faces3 = [f6, f7, f8, f9]
    nodes3 = [n6, n7, n8, n9]
    c3 = TetrahedronCell(faces3, nodes3)

    def test_r_(self):
        self.assertTrue(np.allclose(self.c1.r_, np.array([0.25, 0.25, 0.25])))

    def test_vol(self):
        self.assertEqual(self.c1.vol, 1.0 / 6.0)

    def test_hash_1(self):
        self.assertEqual(hash(self.c1), hash(self.c2))

    def test_hash_2(self):
        self.assertNotEqual(hash(self.c1), hash(self.c3))

if __name__ == '__main__':
    unittest.main()
