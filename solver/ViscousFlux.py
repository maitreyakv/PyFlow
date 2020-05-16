from abc import ABC, abstractmethod

# TODO: Add doc for class
class ViscousFlux(ABC):

    # Constructor for ViscousFlux
    def __init__(self):
        # Superclass constructor
        super().__init__()

    @abstractmethod
    def compute_viscous_flux(self, face, thermo):
        pass
