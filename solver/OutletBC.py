from solver.BC import BC

# TODO: Add doc for class
class OutletBC(BC):

    # Constructor for OutletBC
    def __init__(self, p):
        # Call to superclass constructor
        super().__init__()

        # Save the outlet pressure
        self.p = p

    # Implements the usage of the BC to update a GhostCell
    def apply(self, ghost_cell, interior_cell, boundary_face, thermo, t):
        # Define the reference state using the interior cell
        rho_ref = interior_cell.flow.rho
        c_ref = interior_cell.flow.c

        # Use outflow pressure as the boundary pressure
        boundary_face.flow.p = self.p

        # Compute boundary flow properties
        boundary_face.flow.rho = interior_cell.flow.rho + (boundary_face.flow.p - interior_cell.flow.p) / c_ref**2
        boundary_face.flow.v_ = interior_cell.flow.v_ + boundary_face.n_ * (interior_cell.flow.p - boundary_face.flow.p) / (rho_ref * c_ref)

        # Compute the energy on the boundary
        boundary_face.flow.E = thermo.E(boundary_face.flow)

        # Update conservative variables in the BoundaryFace
        boundary_face.flow.W_[0] = boundary_face.flow.rho
        boundary_face.flow.W_[1:4] = boundary_face.flow.rho * boundary_face.flow.v_
        boundary_face.flow.W_[4] = boundary_face.flow.rho * boundary_face.flow.E

        # Extrapolate conservative variables in the GhostCell from the boundary and interior
        ghost_cell.flow.W_ = 2. * boundary_face.flow.W_ - interior_cell.flow.W_
