import numpy as np

# TODO: Add doc for class
class Node:
    # Constructor for Node
    def __init__(self, x, y, z, id):
        # Save Node coordinates
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

        # Save the ID
        self.id = id

        # Save the position vector of the Node
        self.r_ = np.array([x, y, z], dtype=np.float64)

        # Initialize an empty set of Faces to which the Node belongs to
        self.faces = set()

    # Finds a Face in the set of Faces to which the Nodes belongs to
    # TODO: Test function
    def find_face(self, nodes_face):
        for face_belonged in self.faces:
            if face_belonged.has_nodes(nodes_face):
                return face_belonged
        return None

    # Implement equality function
    def __eq__(self, other):
        # Check of the other Node has the same coordinates
        return all(self.r_ == other.r_)

    # Implement hash function
    def __hash__(self):
        return hash((self.x, self.y, self.z))
