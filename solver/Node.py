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

    # Implement equality function
    def __eq__(self, other):
        # Check if other object is of type Node
        if not isinstance(other, type(self)):
            return False

        # Check of the other Node has the same coordinates
        return np.allclose(self.r_, other.r_)

    # Implement hash function
    def __hash__(self):
        return hash((self.x, self.y, self.z))
