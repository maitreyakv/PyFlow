import numpy as np

from solver.Thermo import Thermo

# TODO: Add doc for class
class CPGThermo(Thermo):

    # Constructor for CPGThermo
    def __init__(self, gamma, Pr, R):
        # Call to superclass constructor
        super().__init__()

        # Save the specific heat capacity ratio
        self.gamma = gamma

        # Save the Prandtl number
        self.Pr = Pr

        # Save the specific gas constant
        self.R = R

    # Pressure EoS
    def p(self, flow):
        return (self.gamma - 1.) * flow.rho * (flow.E - 0.5 * np.linalg.norm(flow.v_)**2)

    # Ideal gas law
    def T(self, flow):
        return flow.p / (flow.rho * self.R)

    # Sutherland's Law
    def mu(self, flow):
        return 1.45 * flow.T**1.5 * 1.e-6 / (flow.T + 110.)
