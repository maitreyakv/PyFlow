import itertools
import numpy as np
from abc import ABC, abstractmethod

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
        self.residual = np.nan * np.ones(5)

    # Adds the convective and viscous fluxes from all the Faces to the residual
    def add_fluxes_to_residual(self):
        for face in self.faces:
            self.residual += (self.Fc_map[face] - self.Fv_map[face]) * face.area

    # Finds the neighbors of the Cell after all Cells have been created
    def find_neighbors(self):
        self.neighbors = [face.other_cell(self) for face in self.faces]

    # Implement hashing function
    def __hash__(self):
        return hash(self.nodes)

    @abstractmethod
    def compute_centroid(self):
        pass

    @abstractmethod
    def compute_volume(self):
        pass
