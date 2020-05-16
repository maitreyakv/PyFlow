import numpy as np

# TODO: Add doc for class
class Flow:

    # Constructor for Flow
    def __init__(self, rho=np.nan, u=np.nan, v=np.nan, w=np.nan, E=np.nan):
        # Initialize vector of conservative variables
        self.W_ = rho * np.array([1., u, v, w, E])

    # Updated Flow variables
    def update(self, thermo):
        # Update density
        self.rho = self.W_[0]

        # Update velocity vector
        self.v_ = self.W_[1:4] / self.rho

        # Update energy
        self.E = self.W_[4] / self.rho

        # Compute pressure
        self.p = thermo.p(self)

        # Compute total enthalpy
        self.H = thermo.H(self)

        # Compute temperature
        self.T = thermo.T(self)

        # Compute dynamic viscosity
        self.mu = thermo.mu(self)

        # Compute thermal conductivity
        self.k = thermo.k(self)

        # Compute speed of sound
        self.c = thermo.c(self)
