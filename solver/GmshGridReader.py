import itertools
import numpy as np
from tqdm import tqdm
import warnings

from solver.Node import Node
from solver.Face import Face
from solver.TetrahedronCell import TetrahedronCell

# TODO: Add doc for class
class GmshGridReader:

    # Constructor for GmshGridReader
    def __init__(self):
        pass

    # TODO: Add doc foc function
    def find_or_create_face(self, nodes_face):
        for node in nodes_face:
            found_face = node.find_face(nodes_face)
            if found_face:
                return found_face
        new_face = Face(*nodes_face)
        self.faces.append(new_face)
        return new_face

    # TODO: Add doc for function
    def read_grid_from_file(self, filename):
        print("reading mesh from Gmsh file...")

        # Initialize a variable for a line in the file
        line = ""

        # Initialize a dictionary for the nodes
        self.nodes = {}

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
            print("reading nodes...")
            for i in tqdm(range(num_nodes)):
                line = fp.readline().strip().split(" ")
                self.nodes[i] = Node(*line[1:], id=i)

            # Read lines until the start of the Elements
            while not "$Elements" in line:
                line = fp.readline()

            # Read the number of Elements
            num_elements = int(fp.readline().strip())

            # Read the Elements and initialize Faces and Cells
            print("reading elements...")
            for i in tqdm(range(num_elements)):
                line = fp.readline().strip().split(" ")

                # Parse Element from file
                if line[1] == '2':
                    face_nodes = (self.nodes[n-1] for n in map(int, line[-3:]))
                    self.faces.append(Face(*face_nodes))
                elif line[1] == '4':
                    nodes_tetra = (self.nodes[n-1] for n in map(int, line[-4:]))
                    faces_tetra = [self.find_or_create_face(nodes_face) for nodes_face in itertools.combinations(nodes_tetra, 3)]
                    self.cells.append(TetrahedronCell(*faces_tetra))

        print("done reading {} nodes, {} faces, {} cells from file".format(len(self.nodes), len(self.faces), len(self.cells)))

        return self.nodes, self.faces, self.cells
