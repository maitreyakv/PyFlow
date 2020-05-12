from solver.BC import BC

# TODO: Add doc for class
class NoSlipAdiabaticWallBC(BC):

    # Constructor for NoSlipAdiabaticWallBC
    def __init__(self):
        # Call to superclass constructor
        super().__init__()

    # Implements the usage of the BC to update a GhostCell
    def apply_bc(self, ghost_cell, interior_cell):
        # Copy the density and energy from the interior Cell to the GhostCell
        ghost_cell.flow.W_[0] = interior_cell.flow.W_[0]
        ghost_cell.flow.W_[4] = interior_cell.flow.W_[4]

        # Reverse the direction of the velocity from the interior Cell to the GhostCell
        ghost_cell.flow.W_[1:4] = -interior_cell.flow.W_[1:4]
