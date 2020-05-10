from solver.Face import Face

# TODO: Add doc for class
class BoundaryFace(Face):
    # Constructor for the BoundaryFace
    def __init__(self, *nodes):
        # Call to superclass constructor
        super().__init__(*nodes)
