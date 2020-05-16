from solver.BC import BC

# TODO: Add doc for class
class InflowBC(BC):

    # Constructor for InflowBC
    def __init__(self, p, rho, v_):
        # Call to superclass constructor
        super().__init__()

        # Save the ambient pressure
        self.p = p

        # Save the ambient density
        self.rho = rho

        # Save the freestream velocity
        self.v_ = v_

    # Implements the usage of the BC to update a GhostCell
    def apply(self, ghost_cell, interior_cell, boundary_face, thermo):
        # Define the reference state using the interior cell
        rho_ref = interior_cell.flow.rho
        c_ref = interior_cell.flow.c

        # Compute the freestream energy
        ghost_cell.flow.p = self.p
        ghost_cell.flow.v_ = self.v_
        ghost_cell.flow.rho = self.rho
        ghost_cell.flow.E = thermo.E(ghost_cell.flow)

        # Set the conservative variables for the GhostCell
        ghost_cell.flow.W_[0] = self.rho
        ghost_cell.flow.W_[1:4] = self.rho * self.v_
        ghost_cell.flow.W_[4] = self.rho * ghost_cell.flow.E

        # Compute boundary pressure
        boundary_face.flow.p = 0.5 * (self.p + interior_cell.flow.p - rho_ref * c_ref * boundary_face.n_.dot(self.v_ - interior_cell.flow.v_))

        # Compute boundary flow properties
        boundary_face.flow.rho = ghost_cell.flow.rho + (boundary_face.flow.p - ghost_cell.flow.p) / c_ref**2
        boundary_face.flow.v_ = ghost_cell.flow.v_ - boundary_face.n_ * (ghost_cell.flow.p - boundary_face.flow.p) / (rho_ref * c_ref)
        boundary_face.flow.E = thermo.E(boundary_face.flow)

        # Set the conservative variables for the boundary
        boundary_face.flow.W_[0] = boundary_face.flow.rho
        boundary_face.flow.W_[1:4] = boundary_face.flow.rho * boundary_face.flow.v_
        boundary_face.flow.W_[4] = boundary_face.flow.rho * boundary_face.flow.E
