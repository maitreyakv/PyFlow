import numpy as np

# TODO: Add doc for class
class Flow:

    # Constructor for Flow
    def __init__(self, rho=np.nan, u=np.nan, v=np.nan, w=np.nan, E=np.nan):
        # Initialize vector of conservative variables
        self.W_ = rho * np.array([1., u, v, w, E])

    # Returns the Flow velocity
    def v_(self):
        return self.W_[1:4] / self.W_[0]
