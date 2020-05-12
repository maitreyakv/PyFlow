from abc import ABC, abstractmethod

# TODO: Add doc for class
class Thermo(ABC):

    # Constructor for Thermo
    def __init__(self):
        pass

    @abstractmethod
    def p(self, flow):
        pass

    @abstractmethod
    def T(self, flow):
        pass

    @abstractmethod
    def mu(self, flow):
        pass
