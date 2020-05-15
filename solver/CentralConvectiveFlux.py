from tqdm import tqdm
import numpy as np

from solver.ConvectiveFlux import ConvectiveFlux

# TODO: Add doc for class
class CentralConvectiveFlux(ConvectiveFlux):

    # Constructor for CentralConvectiveFlux
    def __init__(self, interior_cells, clip=True):
        # Call to superclass constructor
        super().__init__()

        # Create a dictionary mapping Cells to their geometric weights
        self.thetas = {}

        print("initializing central scheme w/ artificial dissipation for convective fluxes...")

        # Compute artificial dissipation geometric quantities for all interior cells
        for cell in tqdm(interior_cells):
            # Compute first order moments
            R = sum([neighbor.r_ - cell.r_ for neighbor in cell.neighbors])

            # Compute second order moments
            I = sum([np.outer(neighbor.r_ - cell.r_, neighbor.r_ - cell.r_) for neighbor in cell.neighbors])

            # Compute matrix of coefficients
            a = np.zeros((3,3))
            for i in range(3):
                for j in range(3):
                    a[i,j] = (-1. ** (i + j)) * np.linalg.det(np.delete(np.delete(I, i, 0), j, 1))

            # Compute determinant of I
            d = np.linalg.det(I)

            # Compute Lagrange multipliers
            lm = (a @ R) / d

            # Compute geometrical weights
            self.thetas[cell] = [1. + lm.dot(neighbor.r_ - cell.r_) for neighbor in cell.neighbors]
            if clip:
                self.thetas[cell] = np.clip(self.thetas[cell], 0., 2.)

    # Implements the computation of the convective flux at a Face
    def compute_convective_flux_at_face(self, face):
        # TODO: Implement function
        pass
