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
        return sum([face.r_.dot(face.normal(self.r_) * face.area) for face in self.faces]) / 3.
