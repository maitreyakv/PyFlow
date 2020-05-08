import itertools
import numpy as np

from solver.Element import Element

# TODO: Add doc for class
class Face(Element):

    # Counter for ID given to each Face instance
    id_iter = itertools.count()

    # Constructor for Face
    def __init__(self, *pts):
        # Initialize the vertices
        self.pts = pts
        if len(pts) < 3:
            raise Exception("Not enough vertices for Face: {}".format(self.id))

        # Assign unique ID to Face using counter
        self.id = next(self.id_iter)

        # Call superclass constructor
        super().__init__()

        # Compute the normal vector
        self.n_ = np.cross(pts[1] - pts[0], pts[2] - pts[0])
        self.n_ /= np.linalg.norm(self.n_)

        # Compute the surface area of the Face
        self.area = abs(0.5 * self.n_.dot(sum([np.cross(pts[i], pts[i+1]) for i in range(len(pts)-1)])))

    # TODO: Add doc for function
    def compute_centroid(self):
        return sum(self.pts) / len(self.pts)

    # Sets the left or right Cell of the Face
    def set_cell(self, cell):
        # TODO: implement and test
        pass

    # Returns the normal unit vector of the Face, pointing away from a reference point
    def normal(self, ref_):
        n_ = self.n_ if self.n_.dot(self.r_ - ref_) > 0.0 else -self.n_
        return n_

    # String representation of Face instance
    def __str__(self):
        # Initialize string representation with ID
        string = "Face {}:\n".format(self.id)

        # Add vertices to string
        string += "\n".join(["  pt.{}: {}".format(i, pt) for i, pt in enumerate(self.pts)])

        # Add normal vector to string
        string += "\n  n_: {}".format(self.n_)

        # Add area to string
        string += "\n  area: {}".format(self.area)

        # Add normal vector to string
        string += "\n  n_: {}".format(self.n_)

        # Add normal vector to string
        string += "\n  c_: {}".format(self.r_)

        # Return string representation
        return string
