from abc import ABC, abstractmethod

# TODO: Add doc for class
class Integrator(ABC):

    # Constructor for Integrator
    def __init__(self):
        # Superclass Constructor
        super().__init__()

    @abstractmethod
    def integrate(self, sim, t, dt=None):
        pass
