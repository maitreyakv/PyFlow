import itertools
import numpy as np
from abc import ABC, abstractmethod

from solver.Element import Element

# TODO: Add doc for class
class Cell(Element, ABC):

    # Constructor for Cell
    def __init__(self, *faces):
        # Initialize the Faces
        self.faces = faces

        # Obtain set of vertices
        self.nodes = tuple({node for face in faces for node in face.nodes})

        # Call superclass constructor
        super().__init__()

        # Compute the volume of the Cell
        self.vol = self.compute_volume()

        # Add the Cell to each Face
        for face in faces:
            face.set_cell(self)

    @abstractmethod
    def compute_centroid(self):
        pass

    @abstractmethod
    def compute_volume(self):
        pass
