from itertools import combinations

from numpy import zeros, fromiter, int64, take, mean, array, concatenate
from tqdm import tqdm

from solver.mesh import node_type, face_type, cell_type
from solver.mesh import create_face, find_or_create_face, create_cell, create_ghost_cell

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
        nodes[i]["rx"] = line[1]
        nodes[i]["ry"] = line[2]
        nodes[i]["rz"] = line[3]

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
            face_node_ids = fromiter(map(int, elem[-3:]), dtype=int64) - 1

            # Create the boundary face
            bc = int(elem[-4])
            next_face_id = create_face(nodes, faces, face_node_ids, next_face_id, bc=bc)

        elif elem[1] == '4': # Create cell
            # Find the IDs of the nodes of the cell
            cell_node_ids = fromiter(map(int, elem[-4:]), dtype=int64) - 1

            cell_face_ids = []
            for face_node_ids in combinations(cell_node_ids, 3):
                found_face_id, next_face_id = find_or_create_face(nodes, faces, array(face_node_ids), next_face_id)
                cell_face_ids.append(found_face_id)

            next_cell_id = create_cell(nodes, faces, cells, cell_node_ids, cell_face_ids, next_cell_id)

    # Close mesh file
    fp.close()

    # Allocate additional space for ghost nodes, ghost_faces, and ghost_cells
    nodes = concatenate((nodes, zeros(num_bndry_faces, dtype=node_type)), axis=0)
    cells = concatenate((cells, zeros(num_bndry_faces, dtype=cell_type)), axis=0)
    faces = concatenate((faces, zeros(3*num_bndry_faces, dtype=face_type)), axis=0)

    # Create ghost cells
    for face in faces:
        if face["bc"] > 0:
            next_face_id, next_cell_id = create_ghost_cell(nodes, faces, cells, face, next_face_id, next_cell_id)

    # Return the nodes, faces, and cells
    return nodes, faces, cells

def write_vtk(filename, nodes, faces, cells, ghost=False):
    fp = open(filename, "w")

    fp.write("# vtk DataFile Version 2.0\ncomment goes here\nASCII\nDATASET UNSTRUCTURED_GRID\n\n")

    fp.write("POINTS {} double\n".format(len(nodes)))
    for node in nodes:
        fp.write("{} {} {}".format(node["rx"], node["ry"], node["rz"]) + "\n")

    if ghost:
        num_cells_to_write = cells.size
    else:
        num_cells_to_write = cells[cells["ghost"] == False].size

    fp.write("\nCELLS {} {}\n".format(num_cells_to_write, 5 * num_cells_to_write))
    for cell in cells:
        if (cell["ghost"] and ghost) or not cell["ghost"]:
            fp.write("4 " + "{} {} {} {}".format(cell["node1"], cell["node2"], cell["node3"], cell["node4"]) + "\n")

    fp.write("\nCELL_TYPES {}\n".format(num_cells_to_write) + "10\n" * num_cells_to_write)

    fp.close()
