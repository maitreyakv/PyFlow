from tqdm import tqdm
from numpy import zeros, delete, clip, outer
from numpy.linalg import det

from solver.ConvectiveFlux import ConvectiveFlux
from solver.BoundaryFace import BoundaryFace

# pressure sensor parameters
# TODO: Maybe not hardcode these, make them options
k2 = 0.5
k4 = 1. / 128.

# TODO: Add doc for class
class CentralConvectiveFlux(ConvectiveFlux):

    # Constructor for CentralConvectiveFlux
    def __init__(self, interior_cells, use_clip=True):
        # Call to superclass constructor
        super().__init__()

        # Create a dictionary mapping Cells to their geometric weights
        self.thetas_map = {}

        print("initializing central scheme w/ artificial dissipation for convective fluxes...")

        # Compute artificial dissipation geometric quantities for all interior cells
        for cell in tqdm(interior_cells):
            # Compute first order moments
            R = sum([neighbor.r_ - cell.r_ for neighbor in cell.neighbors])

            # Compute second order moments
            I = sum([outer(neighbor.r_ - cell.r_, neighbor.r_ - cell.r_) for neighbor in cell.neighbors])

            # Compute matrix of coefficients
            a = zeros((3,3))
            for i in range(3):
                for j in range(3):
                    a[i,j] = (-1. ** (i + j)) * det(delete(delete(I, i, 0), j, 1))

            # Compute determinant of I
            d = det(I)

            # Compute Lagrange multipliers
            lm = (a @ R) / d

            # Compute geometrical weights
            self.thetas_map[cell] = [1. + lm.dot(neighbor.r_ - cell.r_) for neighbor in cell.neighbors]
            if use_clip:
                self.thetas_map[cell] = clip(self.thetas_map[cell], 0., 2.)

        # Initialize a map for the spectral radius of each Cell
        self.spectral_radius_map = {}

        # Initialize a map for the pressure sensor of each Cell
        self.pres_sens_map = {}

        # Initialize a map for the conservative variable Laplacians
        self.lap_map = {}

    # Implements the computation of the convective flux at a Face
    def compute_convective_flux(self, face, thermo):
        # Update conservative variables on interior Face using average of Cell values
        if not isinstance(face, BoundaryFace):
            # TODO: Maybe use interpolation instead of averaging
            face.flow.W_ = 0.5 * (face.left_cell.flow.W_ + face.right_cell.flow.W_)
            face.flow.update(thermo)

        # Compute contravariant velocity at Face
        V = face.flow.v_.dot(face.n_)

        # Compute the convective flux at the Face
        Fc = zeros(5)
        Fc[0] = face.flow.rho * V
        Fc[1:4] = face.flow.rho * face.flow.v_ * V + face.n_ * face.flow.p
        Fc[4] = face.flow.rho * face.flow.H * V

        # Add the convective flux to the map of convective fluxes in the left and right Cell
        face.left_cell.Fc_map[face] = -Fc
        face.right_cell.Fc_map[face] = Fc

    # Updates the spectral radius and pressure sensor maps and conservative variable Laplacians for a Cell
    def prepare_for_artificial_dissipation(self, cell):
        # Compute spectral radius
        self.spectral_radius_map[cell] = sum([(face.flow.v_.dot(face.normal(cell.r_)) + face.flow.c)*face.area for face in cell.faces])

        # Compute pressure sensor
        numer = abs(sum([theta * (neighbor.flow.p - cell.flow.p) for theta, neighbor in zip(self.thetas_map[cell], cell.neighbors)]))
        denom = sum([neighbor.flow.p + cell.flow.p for neighbor in cell.neighbors])
        self.pres_sens_map[cell] = numer / denom

        # Compute conservative variable Laplacians
        self.lap_map[cell] = sum([theta * (neighbor.flow.W_ - cell.flow.W_) for theta, neighbor in zip(self.thetas_map[cell], cell.neighbors)])

    # Computes and add artificial dissipation to a Cell's residual
    def add_artificial_dissipation(self, cell):
        # Initialize dissipation
        D_ = zeros(5)

        # Iterate over Faces and neighbors of Cell
        for face, neighbor, theta in zip(cell.faces, cell.neighbors, self.thetas_map[cell]):
            if not isinstance(face, BoundaryFace):
                # Compute spectral radius at Face
                # TODO: Maybe use interpolation instead of averaging
                spectral_radius = 0.5 * (self.spectral_radius_map[face.left_cell] + self.spectral_radius_map[face.right_cell])

                # Compute pressure sensor coefficients
                eps2 = k2 * max(self.pres_sens_map[cell], self.pres_sens_map[neighbor])
                eps4 = max(0., k4 - eps2)

                # Add contribution to artificial dissipation
                D_ += spectral_radius * eps2 * theta * (neighbor.flow.W_ - cell.flow.W_) - (spectral_radius * eps4 * (self.lap_map[neighbor] - self.lap_map[cell]))

        # Add the artificial dissipation to the Cell's residual
        cell.D_ = D_ # TEMP: Save dissipation for debugging
        cell.Rd += D_
