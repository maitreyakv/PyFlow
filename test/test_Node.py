import unittest
import numpy as np

from solver.Node import Node

class TestNode(unittest.TestCase):

    n1 = Node(1., 2., 3., 101)
    n2 = Node(1., 2., 3., 202)
    n3 = Node(4., 5., 6., 101)

    def test_eq_1(self):
        self.assertTrue(self.n1 == self.n2)

    def test_eq_2(self):
        self.assertFalse(self.n1 == self.n3)

    def test_hash_1(self):
        self.assertTrue(hash(self.n1) == hash(self.n2))

    def test_hash_1(self):
        self.assertFalse(hash(self.n1) == hash(self.n3))

if __name__ == '__main__':
    unittest.main()
