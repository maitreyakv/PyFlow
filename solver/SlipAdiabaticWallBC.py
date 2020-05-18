from solver.BC import BC

# TODO: Add doc for class
class SlipAdiabaticWallBC(BC):

    # Constructor for SlipAdiabaticWallBC
    def __init__(self):
        # Call to superclass constructor
        super().__init__()

    # Implements the usage of the BC to update a GhostCell
    def apply(self, ghost_cell, interior_cell, boundary_face, thermo, t):
        # Copy the density and energy from the interior Cell to the GhostCell
        ghost_cell.flow.W_[0] = interior_cell.flow.W_[0]
        ghost_cell.flow.W_[4] = interior_cell.flow.W_[4]

        # Compute contravariant velocity
        V = interior_cell.flow.v_.dot(boundary_face.n_)

        # Reflect velocity from the interior Cell to the GhostCell
        ghost_cell.flow.v_ = interior_cell.flow.v_ - 2. * V * boundary_face.n_
        ghost_cell.flow.W_[1:4] = ghost_cell.flow.W_[0] * ghost_cell.flow.v_

        # Interpolate conservative variables to the BoundaryFace
        boundary_face.flow.W_ = 0.5 * (ghost_cell.flow.W_ + interior_cell.flow.W_)
