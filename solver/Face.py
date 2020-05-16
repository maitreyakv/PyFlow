import itertools
import numpy as np

from solver.util import fast_cross
from solver.Element import Element

# TODO: Add doc for class
class Face(Element):

    # Constructor for Face
    def __init__(self, *nodes):
        # Initialize the vertices
        self.nodes = tuple(sorted(list(nodes), key=lambda x: hash(x)))

        # Save number of vertices
        self.size = len(nodes)

        # Call superclass constructor
        super().__init__()

        # Compute the normal vector
        self.n_ = fast_cross(nodes[1].r_ - nodes[0].r_, nodes[2].r_ - nodes[0].r_)
        self.n_ /= np.linalg.norm(self.n_)

        # Compute the surface area of the Face
        self.area = 0.5 * np.linalg.norm(fast_cross(nodes[2].r_ - nodes[0].r_, nodes[2].r_ - nodes[1].r_))
        #self.area = abs(0.5 * self.n_.dot(sum([fast_cross(nodes[i].r_, nodes[i+1].r_) for i in range(self.size-1)])))

        # Initialize the left and right cells to None
        self.left_cell = None
        self.right_cell = None

        # Add self to the list of Faces in each Node of the Face
        for node in nodes:
            node.faces.add(self);

    # TODO: Add doc for function
    def compute_centroid(self):
        return sum([node.r_ for node in self.nodes]) / self.size

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

    # Returns whether the Face has these Nodes
    # TODO: Test function
    def has_nodes(self, query_nodes):
        # Return whether all the nodes in the Face
        return all([node in self.nodes for node in query_nodes]) and self.size == len(query_nodes)

    # Returns the other cell
    def other_cell(self, cell):
        return self.left_cell if cell == self.right_cell else self.right_cell

    # Implements equality function
    def __eq__(self, other):
        # Return whether all the nodes are equal
        return self.has_nodes(other.nodes)

    # Implement hash function
    def __hash__(self):
        return hash(self.nodes)
