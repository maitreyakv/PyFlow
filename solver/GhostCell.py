from solver.Cell import Cell

# TODO: Add doc for class
class GhostCell(Cell):

    # Constructor for GhostCell
    def __init__(self, faces, nodes):
        # Call to superclass constructor
        super().__init__(faces, nodes)
