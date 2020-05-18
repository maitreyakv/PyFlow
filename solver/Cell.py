import itertools
from numpy import nan, ones
from abc import ABC, abstractmethod

from solver.util import fast_cross
from solver.Element import Element

# TODO: Add doc for class
class Cell(Element, ABC):

    # Constructor for Cell
    def __init__(self, faces, nodes):
        # Initialize the Faces
        self.faces = faces

        # Obtain set of vertices
        self.nodes = tuple(sorted(list(nodes), key=lambda x: hash(x)))

        # Call superclass constructor
        super().__init__()

        # Compute the volume of the Cell
        self.vol = self.compute_volume()

        # Add the Cell to each Face
        for face in faces:
            face.set_cell(self)

        # Initialize a map from a Face of the Cell to its convective flux
        self.Fc_map = {face: None for face in faces}

        # Initialize a map from a Face of the Cell to its viscous flux
        self.Fv_map = {face: None for face in faces}

        # Initialize a residual for the Cell
        self.Rc = nan * ones(5)
        self.Rd = nan * ones(5)

    # Adds the convective and viscous fluxes from all the Faces to the residual
    def add_fluxes_to_residual(self):
        for face in self.faces:
            self.Rc += self.Fc_map[face] * face.area
            self.Rd += self.Fv_map[face] * face.area

    # Finds the neighbors of the Cell after all Cells have been created
    def find_neighbors(self):
        self.neighbors = [face.other_cell(self) for face in self.faces]

    # Implement hashing function
    def __hash__(self):
        return hash(self.nodes)

    def compute_centroid(self):
        return sum([node.r_ for node in self.nodes]) / 4.0;

    def compute_volume(self):
        # Compute the volume
        return sum([face.r_.dot(face.normal(self.r_) * face.area) for face in self.faces]) / 3.
