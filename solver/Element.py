from abc import ABC, abstractmethod

# TODO: Add doc for class
class Element(ABC):

    # Constructor for Element
    def __init__(self):
        # Compute centroid
        self.r_ = self.compute_centroid()

    @abstractmethod
    def compute_centroid(self):
        pass
