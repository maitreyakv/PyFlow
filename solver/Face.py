import itertools
import numpy as np

from solver.Element import Element

# TODO: Add doc for class
class Face(Element):

    # Counter for ID given to each Face instance
    id_iter = itertools.count()

    # Constructor for Face
    def __init__(self, *pts):
        # Assign unique ID to Face using counter
        self.id = next(self.id_iter)

        # Initialize the vertices
        self.pts = pts
        if len(pts) < 3:
            raise Exception("Not enough vertices for Face: {}".format(self.id))

        # Call superclass constructor
        super().__init__()

        # Compute the normal vector
        self.n_ = np.cross(pts[1] - pts[0], pts[2] - pts[0])
        self.n_ /= np.linalg.norm(self.n_)

        # Compute the surface area of the Face
        self.area = abs(0.5 * self.n_.dot(sum([np.cross(pts[i], pts[i+1]) for i in range(len(pts)-1)])))

        # Initialize the left and righr cells to None
        self.left_cell = None
        self.right_cell = None

    # TODO: Add doc for function
    def compute_centroid(self):
        return sum(self.pts) / len(self.pts)

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

    # Checks if the vertices of the Face are the same as the arguments
    def is_same_vertices(self, *pts):
        # Check if number of vertices are the same
        if len(pts) != len(self.pts):
            return False

        # Check if each query vertex is in the Faces vertices
        for pt in pts:
            if not any([np.allclose(pt, vtx) for vtx in self.pts]):
                return False

        # Return True since query points are vertices of the Face
        return True
