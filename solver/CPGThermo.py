from numpy import sqrt
from numpy.linalg import norm

from solver.Thermo import Thermo

# TODO: Add doc for class
class CPGThermo(Thermo):

    # Constructor for CPGThermo
    def __init__(self, gamma, Pr, R):
        # Call to superclass constructor
        super().__init__()

        # Save the specific heat capacity ratio
        self.gamma_val = gamma

        # Save the Prandtl number
        self.Pr = Pr

        # Save the specific gas constant
        self.R = R

    # Pressure EoS for pressure
    def p(self, flow):
        return (self.gamma_val - 1.) * flow.rho * (flow.E - 0.5 * norm(flow.v_)**2)

    # Pressure EoS for energy
    def E(self, flow):
        return flow.p / ((self.gamma_val - 1.) * flow.rho) + 0.5 * norm(flow.v_)**2

    # Total enthalpy
    def H(self, flow):
        return flow.E + flow.p / flow.rho

    # Ideal gas law for temperature
    def T(self, flow):
        return flow.p / (flow.rho * self.R)

    # Ideal gas law for density
    def rho(self, flow):
        return flow.p / (self.R * flow.T)

    # Speed of sound
    def c(self, flow):
        return sqrt(self.gamma_val * self.R * flow.T)

    # Sutherland's Law
    def mu(self, flow):
        return 1.45 * flow.T**1.5 * 1.e-6 / (flow.T + 110.)

    # Compute thermal conductivity from Prandtl number
    def k(self, flow):
        return self.cp(flow) * flow.mu / self.Pr

    # Return specific heat capacity with constant pressure
    def cp(self, flow):
        return self.gamma_val * self.R / (self.gamma_val - 1.)

    # Return specific heat ratio
    def gamma(self, flow):
        return self.gamma_val
