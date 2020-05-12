import itertools

from solver.Node import Node
from solver.Face import Face
from solver.GhostCell import GhostCell

# TODO: Add doc for class
class BoundaryFace(Face):
    # Constructor for the BoundaryFace
    def __init__(self, *nodes):
        # Call to superclass constructor
        super().__init__(*nodes)

    # Creates a GhostCell
    def create_ghost_cell(self, new_node_id):
        # Get the nodes of the real cell
        nodes_real_cell = self.left_cell.nodes if self.left_cell else self.right_cell.nodes

        # Find the Node in the real Cell not in the BoundaryFace
        node_to_flip = None
        for node in nodes_real_cell:
            if not node in self.nodes:
                node_to_flip = node
                break

        # Create a new Node by reflecting through the BoundaryFace
        new_node_vector = self.reflect_point(node_to_flip.r_)
        new_node = Node(*new_node_vector, new_node_id)

        # Create Faces for GhostCell
        new_faces = [Face(*(node_pair_in_face + (new_node,))) for node_pair_in_face in itertools.combinations(self.nodes, 2)]

        # Create GhostCell
        ghost_cell = GhostCell(tuple(new_faces) + (self,), self.nodes + (new_node,))

        # Assign the GhostCell to the left or right Cell
        if self.left_cell:
            self.right_cell = ghost_cell
        else:
            self.left_cell = ghost_cell

        # Return the Ghost Cell
        return ghost_cell, new_node

    # Helper function to reflect a point about the BoundaryFace
    def reflect_point(self, pt):
        d = -self.n_.dot(self.nodes[0].r_)
        u_ = (pt.dot(self.n_) + d) * self.n_
        v_ = pt - u_
        return -u_ + v_
