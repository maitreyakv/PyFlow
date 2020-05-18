from numpy import int64, float64, dtype, array_equal, sort, min, nonzero

node_type = dtype([("id", int64),
                   ("r_", float64, (3,))])

face_type = dtype([("id",    int64),
                   ("nodes", int64,   (3,)),
                   ("r_",    float64, (3,))])

cell_type = dtype([("id",    int64),
                   ("nodes", int64,   (4,)),
                   ("faces", int64,   (4,)),
                   ("r_",    float64, (3,))])

def create_face(nodes, faces, face_nodes, next_face_id):
    # Assign ID to face
    faces[next_face_id]["id"] = next_face_id

    # Assign the nodes of the face
    faces[next_face_id]["nodes"] = sort(face_nodes)

    # Find the vertices of the face
    vtx_1 = nodes["r_"][face_nodes[0]]
    vtx_2 = nodes["r_"][face_nodes[1]]
    vtx_3 = nodes["r_"][face_nodes[2]]

    # Compute the centroid of the face
    faces[next_face_id]["r_"] = (vtx_1 + vtx_2 + vtx_3) / 3.

    # Increment the number of created faces
    return next_face_id + 1

def find_or_create_face(nodes, faces, face_nodes, next_face_id):
    try:
        return min(nonzero(faces["nodes"] == sort(face_nodes))[0]), next_face_id
    except ValueError:
        next_face_id = create_face(nodes, faces, face_nodes, next_face_id)
        return next_face_id - 1, next_face_id

def create_cell(nodes, faces, cells, cell_nodes, cell_faces, next_cell_id):
    # Assign ID to cell
    cells[next_cell_id]["id"] = next_cell_id

    # Assign the nodes of the cell
    cells[next_cell_id]["nodes"] = sort(cell_nodes)

    # Assign the faces of the cell
    cells[next_cell_id]["faces"] = sort(cell_faces)

    # Find the vertices of the cell
    vtx_1 = nodes["r_"][cell_nodes[0]]
    vtx_2 = nodes["r_"][cell_nodes[1]]
    vtx_3 = nodes["r_"][cell_nodes[2]]
    vtx_4 = nodes["r_"][cell_nodes[3]]

    # Compute the centroid of the cell
    cells[next_cell_id]["r_"] = (vtx_1 + vtx_2 + vtx_3 + vtx_4) / 4.

    # Increment the number of created faces
    return next_cell_id + 1
