import itertools
import numpy as np

from solver.Face import Face
from solver.TetrahedronCell import TetrahedronCell

# TODO: Add doc for class
class GmshGridReader:

    # Constructor for GmshGridReader
    def __init__(self):
        pass

    # TODO: Add doc foc function
    def find_or_create_face(self, pts):
        for face in self.faces:
            if face.is_same_vertices(*pts):
                return face
        self.faces.append(Face(*pts))
        return self.faces[-1]

    # TODO: Add doc for function
    def read_grid_from_file(self, filename):
        # Initialize a variable for a line in the file
        line = ""

        # Initialize a dictionary for the nodes
        nodes = {}

        # Initialze a list for the Faces
        self.faces = []

        # Initialize a list for the Cells
        self.cells = []

        # Attempt to open mesh file
        with open(filename, 'r') as fp:
            # Read lines until the start of the nodes
            while not "$Nodes" in line:
                line = fp.readline()

            # Read the number of nodes in the mesh
            num_nodes = int(fp.readline().strip())

            # Read each node from the file into the dictionary
            for i in range(num_nodes):
                line = fp.readline().strip().split(" ")
                nodes[i+1] = np.array(line[1:], dtype=np.float64)

            # Read lines until the start of the Elements
            while not "$Elements" in line:
                line = fp.readline()

            # Read the number of Elements
            num_elements = int(fp.readline().strip())

            # Read the Elements and initialize Faces and Cells
            for i in range(num_elements):
                line = fp.readline().strip().split(" ")

                # Parse Element from file
                if line[1] == '2':
                    pts = (nodes[n] for n in map(int, line[-3:]))
                    self.faces.append(Face(*pts))
                elif line[1] == '4':
                    pts_tetra = (nodes[n] for n in map(int, line[-4:]))
                    faces_tetra = [self.find_or_create_face(pts) for pts in itertools.combinations(pts_tetra, 3)]
                    self.cells.append(TetrahedronCell(*faces_tetra))
                else:
                    # Unrecognized Element
                    raise Exception("Unrecognized Element in Gmsh file")

        return self.faces, self.cells
