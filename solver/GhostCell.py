from solver.TetrahedronCell import TetrahedronCell

# TODO: Add doc for class
class GhostCell(TetrahedronCell):

    # Constructor for GhostCell
    def __init__(self, faces, nodes):
        # Call to superclass constructor
        super().__init__(faces, nodes)
