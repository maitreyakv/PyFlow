from abc import ABC, abstractmethod

# TODO: Add doc for class
class ConvectiveFlux(ABC):

    # Constructor for ConvectiveFlux
    def __init__(self):
        pass

    @abstractmethod
    def compute_convective_flux(self, face):
        pass
