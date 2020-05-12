from solver.BC import BC

# TODO: Add doc for class
class OutletBC(BC):

    # Constructor for OutletBC
    def __init__(self):
        # Call to superclass constructor
        super().__init__()

    # Implements the usage of the BC to update a GhostCell
    def apply_bc(self, ghost_cell, interior_cell):
        # TEMP: Dummy implementation
        import numpy as np
        ghost_cell.flow.W_[0:5] = np.nan
        # TODO: Implement function
