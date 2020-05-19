from numpy import zeros, float64, where

from solver.thermo import update_flow
from solver.mesh import get_normal, find_face_num_in_cell
from solver.flow import get_velocity

def central_scheme_convective_fluxes(faces, cells, W_faces, W_cells, flow_faces, Fc, thermo="cpg", opts=None):
    # Compute conservative variables on interior faces
    for face in faces:
        # Update conservative variables on interior face using average of cell values
        if face["bc"] == 0:
            face_id = face["id"]
            left_cell_id = face["left_cell"]
            right_cell_id = face["right_cell"]
            W_faces[face_id,:] = 0.5 * ( W_cells[left_cell_id,:] + W_cells[right_cell_id,:] )

    # Update flow variables on all faces
    update_flow(flow_faces, W_faces, thermo=thermo, opts=opts)

    # Compute convective flux across all real faces
    for face in faces:
        if face["bc"] >= 0:
            face_id = face["id"]
            left_cell_id = face["left_cell"]
            right_cell_id = face["right_cell"]
            flow_face = flow_faces[face_id]
            left_cell_face_idx = find_face_num_in_cell(cells[left_cell_id], face) - 1
            right_cell_face_idx = find_face_num_in_cell(cells[right_cell_id], face) - 1

            # Compute contravariant velocity at face
            n_ = get_normal(face)
            v_ = get_velocity(flow_face)
            V = v_.dot(n_)

            # Compute convective flux at face
            flux = zeros(5, dtype=float64)
            flux[0] = flow_face["rho"] * V
            flux[1:4] = flow_face["rho"] * v_ * V + n_ * flow_face["p"]
            flux[4] = flow_face["rho"] * flow_face["H"] * V

            # Add the convective flux to the array of all convective fluxes
            Fc[left_cell_id, left_cell_face_idx, :] = -flux[:]
            Fc[right_cell_id, right_cell_face_idx, :] = flux[:]
