from solver.BC import BC

# TODO: Add doc for class
class NoSlipAdiabaticWallBC(BC):

    # Constructor for NoSlipAdiabaticWallBC
    def __init__(self):
        # Call to superclass constructor
        super().__init__()

    # Implements the usage of the BC to update a GhostCell
    def apply(self, ghost_cell, interior_cell, boundary_face, thermo):
        # Copy the density and energy from the interior Cell to the GhostCell
        ghost_cell.flow.W_[0] = interior_cell.flow.W_[0]
        ghost_cell.flow.W_[4] = interior_cell.flow.W_[4]

        # Reverse the direction of the velocity from the interior Cell to the GhostCell
        ghost_cell.flow.W_[1:4] = -interior_cell.flow.W_[1:4]

        # Set the velocity on the wall to be zero
        boundary_face.flow.W_[1:4] = 0.

        # Copy the density and energy from the interior cell to the wall
        boundary_face.flow.W_[0] =  interior_cell.flow.W_[0]
        boundary_face.flow.W_[4] =  interior_cell.flow.W_[4]
