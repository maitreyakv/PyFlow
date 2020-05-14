from abc import ABC, abstractmethod

from solver.Flow import Flow

# TODO: Add doc for class
class Element(ABC):

    # Constructor for Element
    def __init__(self):
        # Compute centroid
        self.r_ = self.compute_centroid()

        # Declare the Flow field
        self.flow = Flow()

    @abstractmethod
    def compute_centroid(self):
        pass
