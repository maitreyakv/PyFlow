import numpy as np

from solver.Cell import Cell

# TODO: Add doc for class
class TetrahedronCell(Cell):

    # Constructor for the TetrahedronCell
    def __init__(self, *faces):
        # Check that the number of faces is exactly four
        if len(faces) != 4:
            raise Exception("TetrahedronCell needs exactly four faces")

        # Superclass constructor
        super().__init__(*faces)

        # Ensure the TetrahedronCell has exactly four unique vertices
        if len(self.nodes) != 4:
            raise Exception("TetrahedronCell does not have exactly 4 unique vertices")

    def compute_centroid(self):
        return sum([node.r_ for node in self.nodes]) / 4.0;

    def compute_volume(self):
        # Compute the volume
        return abs(np.dot(self.nodes[0].r_ - self.nodes[3].r_, np.cross(self.nodes[1].r_ - self.nodes[3].r_, self.nodes[2].r_ - self.nodes[3].r_))) / 6.0
