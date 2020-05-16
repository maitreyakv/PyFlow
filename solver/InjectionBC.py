from solver.BC import BC

# TODO: Add doc for class
class InjectionBC(BC):

    # Constructor for InjectionBC
    def __init__(self, mass_flow_rate_func, T_inj):
        # Call to superclass constructor
        super().__init__()

        # Save the mass flow rate
        self.mass_flow_rate_func = mass_flow_rate_func

        # Save the injection temperature
        self.T_inj = T_inj

    # Implements the usage of the BC to update a GhostCell
    def apply(self, ghost_cell, interior_cell, boundary_face, thermo, t):
        # Set boundary pressure as interior Cell pressure
        boundary_face.flow.p = interior_cell.flow.p

        # Set boundary temperature as injection temperature
        boundary_face.flow.T = self.T_inj

        # Compute boundary density
        boundary_face.flow.rho = thermo.rho(boundary_face.flow)

        # Compute boundary velocity
        boundary_face.flow.v_ = -boundary_face.n_ * self.mass_flow_rate_func(t) / boundary_face.flow.rho

        # Compute boundary energy
        boundary_face.flow.E = thermo.E(boundary_face.flow)

        # Set the BoundaryFace conservative flow variables
        boundary_face.flow.W_[0] = boundary_face.flow.rho
        boundary_face.flow.W_[1:4] = boundary_face.flow.rho * boundary_face.flow.v_
        boundary_face.flow.W_[4] = boundary_face.flow.rho * boundary_face.flow.E

        # Extrapolate conservative variables to GhostCell
        ghost_cell.flow.W_ = 2. * boundary_face.flow.W_ - interior_cell.flow.W_
