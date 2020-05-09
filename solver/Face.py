import itertools
import numpy as np

from solver.Element import Element

# TODO: Add doc for class
class Face(Element):

    # Constructor for Face
    def __init__(self, *nodes):
        # Initialize the vertices
        self.nodes = tuple(sorted(list(nodes), key=lambda x: hash(x)))
        if len(nodes) < 3:
            raise Exception("Not enough vertices for Face")

        # Call superclass constructor
        super().__init__()

        # Compute the normal vector
        self.n_ = np.cross(nodes[1].r_ - nodes[0].r_, nodes[2].r_ - nodes[0].r_)
        self.n_ /= np.linalg.norm(self.n_)

        # Compute the surface area of the Face
        self.area = abs(0.5 * self.n_.dot(sum([np.cross(nodes[i].r_, nodes[i+1].r_) for i in range(len(nodes)-1)])))

        # Initialize the left and righr cells to None
        self.left_cell = None
        self.right_cell = None

    # TODO: Add doc for function
    def compute_centroid(self):
        return sum([node.r_ for node in self.nodes]) / len(self.nodes)

    # Sets the left or right Cell of the Face
    def set_cell(self, cell):
        # If the Cell is on the right side of the face, set is as the right Cell
        if self.n_.dot(cell.r_ - self.r_) > 0.0:
            self.right_cell = cell
        # Otherwise set the Cell as the left Cell
        else:
            self.left_cell = cell

    # Returns the normal unit vector of the Face, pointing away from a reference point
    def normal(self, ref_):
        n_ = self.n_ if self.n_.dot(self.r_ - ref_) > 0.0 else -self.n_
        return n_

    # Implements equality function
    def __eq__(self, other):
        # Check if other object is of type Face
        if not isinstance(other, type(self)):
            return False

        # Check if the other Face has the same number of vertices
        if len(self.nodes) != len(other.nodes):
            return False

        # Return whether all the nodes are equal
        return all([node in other.nodes for node in self.nodes])

    # Implement hash function
    def __hash__(self):
        return hash(self.nodes)
