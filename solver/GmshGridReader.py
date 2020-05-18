from itertools import combinations
from tqdm import tqdm

from solver.Node import Node
from solver.Face import Face
from solver.BoundaryFace import BoundaryFace
from solver.TetrahedronCell import TetrahedronCell
from solver.InjectionBC import InjectionBC
from solver.OutletBC import OutletBC
from solver.SlipAdiabaticWallBC import SlipAdiabaticWallBC
from solver.NoSlipAdiabaticWallBC import NoSlipAdiabaticWallBC

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

        # Dictionary for the groups and their IDs
        groups = {}

        # Attempt to open mesh file
        with open(filename, 'r') as fp:
            # Read until the start of the groups
            while not "$PhysicalNames" in line:
                line = fp.readline()

            # Read the number of groups in the mesh
            num_groups = int(fp.readline().strip())

            # Read each group from the file into the dictionary
            print("reading groups...")
            for i in range(num_groups):
                line = fp.readline().strip().split(" ")
                if line[0] == "2":
                    groups[int(line[1])] = line[2]

            # Read lines until the start of the nodes
            while not "$Nodes" in line:
                line = fp.readline()

            # Read the number of nodes in the mesh
            num_nodes = int(fp.readline().strip())

            # Read each node from the file into the dictionary
            print("reading nodes...")
            for i in range(num_nodes):
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
                    group = groups[int(line[-4])]
                    if "injection" in group:
                        # TEMP: Hardcoded inlet variables, needs proper IO
                        bc = InjectionBC(lambda t: 1. * 101325. / (287.058 * 273.), 273.)
                    elif "outlet" in group:
                        # TEMP: Hardcoded outlet pressure, needs proper IO
                        bc = OutletBC(101325.) # * 0.5)
                    elif "wall" in group:
                        bc = SlipAdiabaticWallBC()
                    else:
                        raise Exception("Unrecognized BC in Gmsh file")
                    self.faces.append(BoundaryFace(*face_nodes, bc=bc))
                elif line[1] == '4':
                    nodes_tetra = [self.nodes[n-1] for n in map(int, line[-4:])]
                    faces_tetra = [self.find_or_create_face(nodes_face) for nodes_face in combinations(nodes_tetra, 3)]
                    self.cells.append(TetrahedronCell(faces_tetra, nodes_tetra))

        print("done reading {} nodes, {} faces, {} cells from file".format(len(self.nodes), len(self.faces), len(self.cells)))

        print("creating ghost cells...")

        # Create GhostCells
        for face in tqdm([face for face in self.faces if isinstance(face, BoundaryFace)]):
            new_node_id = len(self.nodes)
            ghost_cell, new_node = face.create_ghost_cell(new_node_id)
            self.cells.append(ghost_cell)
            self.nodes[new_node_id] = new_node

        # Determine all neighbors for each cell
        for cell in self.cells:
            cell.find_neighbors()

        return self.nodes, self.faces, self.cells
