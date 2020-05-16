import numpy as np

from solver.Gradient import Gradient

# TODO: Add doc for class
class GreenGaussGradient(Gradient):

    # Constructor for GreenGaussGradient
    def __init__(self):
        # Superclass constructor
        super().__init__()

    # Implements computation of the velocity and temperature gradients
    def compute_gradients(self, cell):
        # Initialize Cell gradients to zero
        cell.flow.grad_u = np.zeros(3)
        cell.flow.grad_v = np.zeros(3)
        cell.flow.grad_w = np.zeros(3)
        cell.flow.grad_T = np.zeros(3)

        # Compute gradients by looping over the Faces of the Cell
        for face in cell.faces:
            cell.flow.grad_u += face.flow.v_[0] * face.normal(cell.r_) * face.area
            cell.flow.grad_v += face.flow.v_[1] * face.normal(cell.r_) * face.area
            cell.flow.grad_w += face.flow.v_[2] * face.normal(cell.r_) * face.area
            cell.flow.grad_T += face.flow.T * face.normal(cell.r_) * face.area

        # Normalize by the cell volume to finish gradient computation
        cell.flow.grad_u /= cell.vol
        cell.flow.grad_v /= cell.vol
        cell.flow.grad_w /= cell.vol
        cell.flow.grad_T /= cell.vol
