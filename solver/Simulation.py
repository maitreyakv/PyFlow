from numpy import zeros, sqrt

from solver.GhostCell import GhostCell
from solver.BoundaryFace import BoundaryFace

# TODO: Add doc for class
class Simulation:

    # Constructor for Simulation
    def __init__(self, faces, cells):
        self.faces = faces
        self.cells = cells
        self.interior_cells = [cell for cell in cells if not isinstance(cell, GhostCell)]
        self.ghost_cells = [cell for cell in cells if isinstance(cell, GhostCell)]
        self.boundary_faces = [face for face in faces if isinstance(face, BoundaryFace)]

        # Initialize Thermo
        from solver.CPGThermo import CPGThermo
        self.thermo = CPGThermo(1.4, 0.72, 287.058)

        # Initialize Convective Flux Scheme
        from solver.CentralConvectiveFlux import CentralConvectiveFlux
        self.convective_flux = CentralConvectiveFlux(self.interior_cells)

        # Initialize Gradient computer
        from solver.LeastSquaresGradient import LeastSquaresGradient
        self.gradient = LeastSquaresGradient(self.interior_cells)

        # Initialize Viscous Flux Scheme
        from solver.GradientAvgViscousFlux import GradientAvgViscousFlux
        self.viscous_flux = GradientAvgViscousFlux()

    # TODO: Add doc for function
    def compute_residuals(self, t):
        # Update all variables in interior cells and set residuals to zero
        for interior_cell in self.interior_cells:
            interior_cell.flow.update(self.thermo)
            interior_cell.Rc = zeros(5)
            interior_cell.Rd = zeros(5)

        # Apply BCs to GhostCells and BoundaryFaces and update all variables on BoundaryFaces
        for boundary_face in self.boundary_faces:
            boundary_face.apply_bc(self.thermo, t)
            boundary_face.flow.update(self.thermo)

        # Update all variables in GhostCells
        for ghost_cell in self.ghost_cells:
            ghost_cell.flow.update(self.thermo)

        # Compute the convective fluxes
        #print("computing convective fluxes...")
        for face in self.faces:
            self.convective_flux.compute_convective_flux(face, self.thermo)

        # Update Cell's spectral radius of convective flux jacobian
        for cell in self.interior_cells:
            self.convective_flux.prepare_for_artificial_dissipation(cell)

        # Add the artificial dissipation
        for cell in self.interior_cells:
            self.convective_flux.add_artificial_dissipation(cell)

        # Compute velocity and temperature gradients
        for cell in self.interior_cells:
            self.gradient.compute_gradients(cell)

        # Compute the viscous fluxes
        #print("computing viscous fluxes...")
        for face in self.faces:
            self.viscous_flux.compute_viscous_flux(face, self.thermo)

        # Add the fluxes to the residuals and compute L2 norm of residuals
        residual_L2_norm = zeros(5)
        for cell in self.interior_cells:
            cell.add_fluxes_to_residual()
            residual_L2_norm += (cell.Rc - cell.Rd)**2
        residual_L2_norm = sqrt(residual_L2_norm)
        print("residual L2 norm = {}".format(residual_L2_norm))

    # TODO: Add doc for function
    def integrate_step(self, dt):
        # TEMP: Perform Euler integration
        for cell in self.interior_cells:
            cell.flow.W_ += -dt * (cell.Rc - cell.Rd) / cell.vol

    # TODO: Add doc for function
    def prepare_for_save(self):
        # Update all variables in all cells before writing to file
        for cell in self.cells:
            cell.flow.update(self.thermo)
