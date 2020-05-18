from numpy import vstack, ones, array
from numpy.linalg import lstsq

from solver.Gradient import Gradient

# TODO: Add doc for class
class LeastSquaresGradient(Gradient):

    # Constructor for LeastSquaresGradient
    def __init__(self, cells):
        # Superclass constructor
        super().__init__()

        # Create a map for the edge matrix
        self.edge_matrix_map = {}

        # Create a map for the weight matrix
        self.weight_map = {}

        # Compute the edge matrix and weight for each cell
        for cell in cells:
            self.edge_matrix_map[cell] = vstack((neighbor.r_ - cell.r_ for neighbor in cell.neighbors))
            # TODO: Look into alternative weights
            self.weight_map[cell] = ones(len(cell.faces))

    # Implements computation of the velocity and temperature gradients
    def compute_gradients(self, cell):
        # Get weight map for Cell
        weight_map = self.weight_map[cell]

        # Assemble difference vectors
        diff_u = weight_map * array([neighbor.flow.v_[0] - cell.flow.v_[0] for neighbor in cell.neighbors])
        diff_v = weight_map * array([neighbor.flow.v_[1] - cell.flow.v_[1] for neighbor in cell.neighbors])
        diff_w = weight_map * array([neighbor.flow.v_[2] - cell.flow.v_[2] for neighbor in cell.neighbors])
        diff_T = weight_map * array([neighbor.flow.T - cell.flow.T for neighbor in cell.neighbors])

        # Get edge matrix for Cell
        edge_matrix = self.edge_matrix_map[cell]

        # Compute derivatives
        cell.flow.grad_u_, _, _, _ = lstsq(edge_matrix, diff_u, rcond=None)
        cell.flow.grad_v_, _, _, _ = lstsq(edge_matrix, diff_v, rcond=None)
        cell.flow.grad_w_, _, _, _ = lstsq(edge_matrix, diff_w, rcond=None)
        cell.flow.grad_T_, _, _, _ = lstsq(edge_matrix, diff_T, rcond=None)
