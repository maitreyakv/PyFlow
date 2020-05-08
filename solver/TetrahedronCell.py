import numpy as np

from solver.Cell import Cell

# TODO: Add doc for class
class TetrahedronCell(Cell):

    # Constructor for the TetrahedronCell
    def __init__(self, *faces):
        # Check that the number of faces is exactly four
        if len(faces) != 4:
            raise Exception("TetrahedronCell {} needs exactly four faces".format(self.id))

        # Superclass constructor
        super().__init__(*faces)

        # Ensure the TetrahedronCell has exactly four unique vertices
        if len(self.pts) != 4:
            raise Exception("TetrahedronCell {} does not have exactly 4 unique vertices".format(self.id))

    def compute_centroid(self):
        return sum(self.pts) / 4.0;

    def compute_volume(self):
        # Compute the volume
        return abs(np.dot(self.pts[0] - self.pts[3], np.cross(self.pts[1] - self.pts[3], self.pts[2] - self.pts[3]))) / 6.0
