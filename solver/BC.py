from abc import ABC, abstractmethod

# TODO: Add doc for class
class BC(ABC):

    # Constructor for BC
    def __init__(self):
        pass

    @abstractmethod
    def apply_bc(self, ghost_cell, interior_cell):
        pass
