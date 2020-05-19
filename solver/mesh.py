from numpy import int64, float64, dtype, array_equal, sort, min, nonzero, array, amax, any
from numpy.linalg import norm
#from numba import jit

from solver.util import fast_cross

node_type = dtype([("id", int64),
                   ("rx", float64),
                   ("ry", float64),
                   ("rz", float64)])

face_type = dtype([("id",         int64),
                   ("node1",      int64),
                   ("node2",      int64),
                   ("node3",      int64),
                   ("rx",         float64),
                   ("ry",         float64),
                   ("rz",         float64),
                   ("nx",         float64),
                   ("ny",         float64),
                   ("nz",         float64),
                   ("area",       float64),
                   ("left_cell",  int64),
                   ("right_cell", int64),
                   ("bc",         int64)])

cell_type = dtype([("id",    int64),
                   ("node1", int64),
                   ("node2", int64),
                   ("node3", int64),
                   ("node4", int64),
                   ("face1", int64),
                   ("face2", int64),
                   ("face3", int64),
                   ("face4", int64),
                   ("rx",    float64),
                   ("ry",    float64),
                   ("rz",    float64),
                   ("vol",   float64),
                   ("ghost", bool)])

#@jit(nopython=True)
def get_normal(face):
    return array( [ face["nx"], face["ny"], face["nz"] ] )

#@jit(nopython=True)
def get_centroid(elem):
    return array( [ elem["rx"], elem["ry"], elem["rz"] ] )

#@jit(nopython=True)
def compute_area_and_normal(nodes, node1, node2, node3):
    vtx_1 = get_centroid(node1)
    vtx_2 = get_centroid(node2)
    vtx_3 = get_centroid(node3)
    area_ = fast_cross(vtx_2 - vtx_1, vtx_3 - vtx_1)
    n_ = area_ / norm(area_)
    return norm(area_), n_[0], n_[1], n_[2]

#@jit(nopython=True)
def compute_volume(faces, face1, face2, face3, face4):
    r1_ = get_centroid(face1)
    r2_ = get_centroid(face2)
    r3_ = get_centroid(face3)
    r4_ = get_centroid(face4)
    n1_ = get_normal(face1)
    n2_ = get_normal(face2)
    n3_ = get_normal(face3)
    n4_ = get_normal(face4)
    area1 = face1["area"]
    area2 = face2["area"]
    area3 = face3["area"]
    area4 = face4["area"]
    return ( r1_.dot(n1_) * area1
           + r2_.dot(n2_) * area2
           + r3_.dot(n3_) * area3
           + r4_.dot(n4_) * area4 ) / 3.

def create_face(nodes, faces, face_node_ids, next_face_id, bc=0):
    # Assign ID to face
    faces[next_face_id]["id"] = next_face_id

    # Assign the nodes of the face
    face_node_ids = sort(face_node_ids)
    node1 = nodes[face_node_ids[0]]
    node2 = nodes[face_node_ids[1]]
    node3 = nodes[face_node_ids[2]]
    faces[next_face_id]["node1"] = node1["id"]
    faces[next_face_id]["node2"] = node2["id"]
    faces[next_face_id]["node3"] = node3["id"]

    # Compute the centroid of the face
    faces[next_face_id]["rx"] = ( node1["rx"] + node2["rx"] + node3["rx"] ) / 3.
    faces[next_face_id]["ry"] = ( node1["ry"] + node2["ry"] + node3["ry"] ) / 3.
    faces[next_face_id]["rz"] = ( node1["rz"] + node2["rz"] + node3["rz"] ) / 3.

    # Compute area and normal vector
    area, nx, ny, nz = compute_area_and_normal(nodes, node1, node2, node3)
    faces[next_face_id]["area"] = area
    faces[next_face_id]["nx"] = nx
    faces[next_face_id]["ny"] = ny
    faces[next_face_id]["nz"] = nz

    # Set the boundary condition
    faces[next_face_id]["bc"] = bc

    # Increment the number of created faces
    return next_face_id + 1

def find_or_create_face(nodes, faces, face_node_ids, next_face_id):
    face_node_ids = sort(face_node_ids)
    node1 = nodes[face_node_ids[0]]
    node2 = nodes[face_node_ids[1]]
    node3 = nodes[face_node_ids[2]]

    mask = (faces["node1"] == node1["id"]) & (faces["node2"] == node2["id"]) & (faces["node3"] == node3["id"])
    if any(mask):
        return min(nonzero(mask)[0]), next_face_id
    else:
        next_face_id = create_face(nodes, faces, face_node_ids, next_face_id)
        return next_face_id - 1, next_face_id

def create_cell(nodes, faces, cells, cell_node_ids, cell_face_ids, next_cell_id, ghost=False):
    # Assign ID to cell
    cells[next_cell_id]["id"] = next_cell_id

    # Assign the nodes of the cell
    cell_node_ids = sort(cell_node_ids)
    node1 = nodes[cell_node_ids[0]]
    node2 = nodes[cell_node_ids[1]]
    node3 = nodes[cell_node_ids[2]]
    node4 = nodes[cell_node_ids[3]]
    cells[next_cell_id]["node1"] = node1["id"]
    cells[next_cell_id]["node2"] = node2["id"]
    cells[next_cell_id]["node3"] = node3["id"]
    cells[next_cell_id]["node4"] = node4["id"]

    # Assign the faces of the cell
    cell_face_ids = sort(cell_face_ids)
    face1 = faces[cell_face_ids[0]]
    face2 = faces[cell_face_ids[1]]
    face3 = faces[cell_face_ids[2]]
    face4 = faces[cell_face_ids[3]]
    cells[next_cell_id]["face1"] = face1["id"]
    cells[next_cell_id]["face2"] = face2["id"]
    cells[next_cell_id]["face3"] = face3["id"]
    cells[next_cell_id]["face4"] = face4["id"]

    # Compute cell volume
    cells[next_cell_id]["vol"] = compute_volume(faces, face1, face2, face3, face4)

    # Compute the centroid of the cell
    cells[next_cell_id]["rx"] = ( node1["rx"] + node2["rx"] + node3["rx"] + node4["rx"] ) / 4.
    cells[next_cell_id]["ry"] = ( node1["ry"] + node2["ry"] + node3["ry"] + node4["ry"] ) / 4.
    cells[next_cell_id]["rz"] = ( node1["rz"] + node2["rz"] + node3["rz"] + node4["rz"] ) / 4.

    # Sets the ghost status of the cell
    cells[next_cell_id]["ghost"] = ghost

    # Associate the cell with its faces
    set_cell_of_face(face1, cells[next_cell_id])
    set_cell_of_face(face2, cells[next_cell_id])
    set_cell_of_face(face3, cells[next_cell_id])
    set_cell_of_face(face4, cells[next_cell_id])

    # Increment the number of created faces
    return next_cell_id + 1

def set_cell_of_face(face, cell):
    n_      = get_normal(face)
    r_face_ = get_centroid(face)
    r_cell_ = get_centroid(cell)
    if not face["bc"]:
        if n_.dot(r_cell_ - r_face_) > 0.:
            face["right_cell"] = cell["id"]
        else:
            face["left_cell"] = cell["id"]
    else:
        if cell["ghost"]:
            return
        if n_.dot(r_cell_ - r_face_) > 0.:
            n_ = -n_
            face["nx"] = -face["nx"]
            face["ny"] = -face["ny"]
            face["nz"] = -face["nz"]
        face["left_cell"] = cell["id"]

def reflect_point(v1_, face, nodes):
    n_ = get_normal(face)
    node1 = nodes[face["node1"]]
    r_node1_ = get_centroid(node1)
    d = -n_.dot(r_node1_)
    u_ = (v1_.dot(n_) + d) * n_
    v_ = v1_ - u_
    return -u_ + v_

def create_ghost_cell(nodes, faces, cells, face, next_face_id, next_cell_id):
    real_node_1_id = cells[face["left_cell"]]["node1"]
    real_node_2_id = cells[face["left_cell"]]["node2"]
    real_node_3_id = cells[face["left_cell"]]["node3"]
    real_node_4_id = cells[face["left_cell"]]["node4"]
    real_node_ids = [real_node_1_id, real_node_2_id, real_node_3_id, real_node_4_id]
    node_to_flip_id = None
    for node_id in real_node_ids:
        if not (node_id == face["node1"] or node_id == face["node2"] or node_id == face["node3"]):
            node_to_flip_id = node_id
    node_to_flip = nodes[node_to_flip_id]
    point_to_flip = get_centroid(node_to_flip)
    new_point = reflect_point(point_to_flip, face, nodes)
    ghost_node_id = amax(nodes["id"]) + 1
    nodes[ghost_node_id]["id"] = ghost_node_id
    nodes[ghost_node_id]["rx"] = new_point[0]
    nodes[ghost_node_id]["ry"] = new_point[1]
    nodes[ghost_node_id]["rz"] = new_point[2]
    ghost_face_1_node_ids = [ face["node1"], face["node2"], ghost_node_id ]
    ghost_face_2_node_ids = [ face["node2"], face["node3"], ghost_node_id ]
    ghost_face_3_node_ids = [ face["node1"], face["node3"], ghost_node_id ]
    next_face_id = create_face(nodes, faces, ghost_face_1_node_ids, next_face_id, bc=-1)
    next_face_id = create_face(nodes, faces, ghost_face_2_node_ids, next_face_id, bc=-1)
    next_face_id = create_face(nodes, faces, ghost_face_3_node_ids, next_face_id, bc=-1)
    ghost_cell_node_ids = [ face["node1"], face["node2"], face["node3"], ghost_node_id ]
    ghost_cell_face_ids = [ next_face_id - 1, next_face_id - 2, next_face_id - 3, face["id"] ]
    next_cell_id = create_cell(nodes, faces, cells, ghost_cell_node_ids, ghost_cell_face_ids, next_cell_id, ghost=True)
    face["right_cell"] = next_cell_id - 1
    return next_face_id, next_cell_id
