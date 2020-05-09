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

    f1 = Face(n1, n2, n3)
    f2 = Face(n1, n2, n4)
    f3 = Face(n1, n3, n4)
    f4 = Face(n2, n3, n4)
    faces = [f1, f2, f3, f4]
    nodes = [n1, n2, n3, n4]
    c1 = TetrahedronCell(faces, nodes)

    def test_r_(self):
        self.assertTrue(np.allclose(self.c1.r_, np.array([0.25, 0.25, 0.25])))

    def test_vol(self):
        self.assertEqual(self.c1.vol, 1.0 / 6.0)

if __name__ == '__main__':
    unittest.main()
