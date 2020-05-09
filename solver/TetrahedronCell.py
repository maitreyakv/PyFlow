import numpy as np

from solver.util import fast_cross
from solver.Cell import Cell

# TODO: Add doc for class
class TetrahedronCell(Cell):

    # Constructor for the TetrahedronCell
    def __init__(self, *faces):
        # Superclass constructor
        super().__init__(*faces)

    def compute_centroid(self):
        return sum([node.r_ for node in self.nodes]) / 4.0;

    def compute_volume(self):
        # Compute the volume
        a = self.nodes[0].r_
        b = self.nodes[1].r_
        c = self.nodes[2].r_
        d = self.nodes[3].r_
        return abs(np.dot(a - d, fast_cross(b - d, c - d))) / 6.0
