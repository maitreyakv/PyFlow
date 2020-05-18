from itertools import combinations

from numpy import zeros, fromiter, int64, take, mean, array
from tqdm import tqdm

from solver.mesh import node_type, face_type, cell_type
from solver.mesh import create_face, find_or_create_face, create_cell

# TODO: Cleanup and add doc

def read_until(fp, stop_str):
    """Reads lines from a file pointer until the line contains a stop string"""
    line = fp.readline()
    while not stop_str in line:
        line = fp.readline()

def read_gmsh_2(filename):
    """Reads an unstructured mesh from a Gmsh 2.0 file"""
    # Open mesh file
    fp = open(filename, 'r')

    # Read until the start of the groups
    read_until(fp, "$PhysicalNames")

    # Read the number of groups in the mesh
    num_groups = int( fp.readline().strip() )

    # Read each group from the file into the dictionary
    groups = {}
    for i in range(num_groups):
        line = fp.readline().strip().split(" ")
        if line[0] == "2":
            groups[ int(line[1]) ] = line[2]

    # Read lines until the start of the nodes
    read_until(fp, "$Nodes")

    # Read the number of nodes in the mesh
    num_nodes = int( fp.readline().strip() )

    # Allocate array for the Nodes
    nodes = zeros(num_nodes, dtype=node_type)

    # Read each node from the file into the dictionary
    for i in range(num_nodes):
        line = fp.readline().strip().split(" ")
        nodes[i]["id"] = i
        nodes[i]["r_"] = line[1:]

    # Read lines until the start of the Elements
    read_until(fp, "$Elements")

    # Read the number of Elements
    num_elements = int( fp.readline().strip() )

    # Read all the lines for the elements into a list
    elem_lines = []
    for _ in range(num_elements):
        elem_lines.append( fp.readline().strip().split(" ") )

    # Determine number of boudary faces
    num_bndry_faces = sum(line[1] == '2' for line in elem_lines)

    # Determine number of cells
    num_cells = sum(line[1] == '4' for line in elem_lines)

    # Determine total number of faces (boundary and internal)
    num_faces = num_bndry_faces + int((4 * num_cells - num_bndry_faces) / 2)

    # Allocate arrays for the faces and cells
    faces = zeros(num_faces, dtype=face_type)
    cells = zeros(num_cells, dtype=cell_type)

    # Initialize counters for the number of created faces and cells
    next_face_id = 0
    next_cell_id = 0

    # Parse the elements and create faces and cells
    for elem in tqdm(elem_lines):
        if elem[1] == '2': # Create boundary face
            # Find the IDs of the nodes of the boundary face
            face_nodes = fromiter(map(int, elem[-3:]), dtype=int64) - 1

            # Create the boundary face
            next_face_id = create_face(nodes, faces, face_nodes, next_face_id)

        elif elem[1] == '4': # Create cell
            # Find the IDs of the nodes of the cell
            cell_nodes = fromiter(map(int, elem[-4:]), dtype=int64) - 1

            cell_faces = []
            for face_nodes in combinations(cell_nodes, 3):
                found_face, next_face_id = find_or_create_face(nodes, faces, array(face_nodes), next_face_id)
                cell_faces.append(found_face)

            next_cell_id = create_cell(nodes, faces, cells, cell_nodes, cell_faces, next_cell_id)

    # Close mesh file
    fp.close()

    # Return the nodes, faces, and cells
    return nodes, faces, cells

def write_vtk(filename, nodes, faces, cells):
    fp = open(filename, "w")

    fp.write("# vtk DataFile Version 2.0\ncomment goes here\nASCII\nDATASET UNSTRUCTURED_GRID\n\n")

    fp.write("POINTS {} double\n".format(len(nodes)))
    for id in range(nodes.size):
        fp.write(" ".join(map(str, nodes[id]["r_"].tolist())) + "\n")

    fp.write("\nCELLS {} {}\n".format(cells.size, 5 * cells.size))
    for id in range(cells.size):
        fp.write("4 " + " ".join(map(str, cells[id]["nodes"])) + "\n")

    fp.write("\nCELL_TYPES {}\n".format(cells.size) + "10\n" * cells.size)

    fp.close()
