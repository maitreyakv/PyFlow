from abc import ABC, abstractmethod

# TODO: Add doc for class
class Gradient(ABC):

    # Constructor for Gradient
    def __init__(self):
        # Superclass constructor
        super().__init__()

    @abstractmethod
    def compute_gradients(self, cell):
        pass
