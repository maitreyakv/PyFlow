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
    def E(self, flow):
        pass

    @abstractmethod
    def H(self, flow):
        pass

    @abstractmethod
    def T(self, flow):
        pass

    @abstractmethod
    def rho(self, flow):
        pass

    @abstractmethod
    def c(self, flow):
        pass

    @abstractmethod
    def mu(self, flow):
        pass

    @abstractmethod
    def k(self, flow):
        pass

    @abstractmethod
    def cp(self, flow):
        pass

    @abstractmethod
    def gamma(self, flow):
        pass
