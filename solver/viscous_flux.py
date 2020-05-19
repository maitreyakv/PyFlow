from numpy import zeros, float64
from numpy.linalg import norm

from solver.mesh import get_centroid, get_normal, find_face_num_in_cell
from solver.flow import get_velocity

# TODO: Add doc and cleanup

def gradient_avg_viscous_fluxes(faces, cells, flow_faces, flow_cells, grads, Fv):
    # Compute viscous fluxes across all real faces
    for face in faces:
        if face["bc"] >= 0:
            face_id = face["id"]
            left_cell_id = face["left_cell"]
            right_cell_id = face["right_cell"]
            flow_face = flow_faces[face_id]
            left_cell_face_idx = find_face_num_in_cell(cells[left_cell_id], face) - 1
            right_cell_face_idx = find_face_num_in_cell(cells[right_cell_id], face) - 1

            # If boundary face, then use interior cell gradients as face gradients
            if face["bc"] > 0:
                inter_cell_id = face["left_cell"]
                grad_u_ = grads[inter_cell_id,0,:]
                grad_v_ = grads[inter_cell_id,1,:]
                grad_w_ = grads[inter_cell_id,2,:]
                grad_T_ = grads[inter_cell_id,3,:]
            else: # Compute gradients at faces from left and right cell gradients
                # Compute simple average of left and right cell gradients
                grad_u_avg_ = 0.5 * ( grads[left_cell_id,0,:] + grads[right_cell_id,0,:] )
                grad_v_avg_ = 0.5 * ( grads[left_cell_id,1,:] + grads[right_cell_id,1,:] )
                grad_w_avg_ = 0.5 * ( grads[left_cell_id,2,:] + grads[right_cell_id,2,:] )
                grad_T_avg_ = 0.5 * ( grads[left_cell_id,3,:] + grads[right_cell_id,3,:] )

                # Compute distance betwen left and right cell centroids
                r_left_cell_ = get_centroid(cells[left_cell_id])
                r_right_cell_ = get_centroid(cells[right_cell_id])
                l = norm(r_right_cell_ - r_left_cell_)

                # Compute unit vector along this line
                t_ = (r_right_cell_ - r_left_cell_) / l

                # Get left and right cell flow states
                flow_left_cell = flow_cells[left_cell_id]
                flow_right_cell = flow_cells[right_cell_id]

                # Compute direcitonal derivatives along line connecting cell centroids
                dudl = (flow_right_cell["u"] - flow_left_cell["u"]) / l
                dvdl = (flow_right_cell["v"] - flow_left_cell["v"]) / l
                dwdl = (flow_right_cell["w"] - flow_left_cell["w"]) / l
                dTdl = (flow_right_cell["T"] - flow_left_cell["T"]) / l

                # Compute modified average of gradients
                grad_u_ = grad_u_avg_ - (grad_u_avg_.dot(t_) - dudl) * t_
                grad_v_ = grad_v_avg_ - (grad_v_avg_.dot(t_) - dvdl) * t_
                grad_w_ = grad_w_avg_ - (grad_w_avg_.dot(t_) - dwdl) * t_
                grad_T_ = grad_T_avg_ - (grad_T_avg_.dot(t_) - dTdl) * t_

            # Compute divergence of velocity
            div_ = grad_u_[0] + grad_v_[1] + grad_w_[2]

            # Assemble stress tensor
            tau = zeros((3,3), dtype=float64)
            tau[0,0] = 2. * flow_face["mu"] * (grad_u_[0] - div_ / 3.)
            tau[1,1] = 2. * flow_face["mu"] * (grad_v_[1] - div_ / 3.)
            tau[2,2] = 2. * flow_face["mu"] * (grad_w_[2] - div_ / 3.)
            tau[0,1] = tau[1,0] = flow_face["mu"] * (grad_u_[1] + grad_v_[0])
            tau[0,2] = tau[2,0] = flow_face["mu"] * (grad_u_[2] + grad_w_[0])
            tau[1,2] = tau[2,1] = flow_face["mu"] * (grad_v_[2] + grad_w_[1])

            # Compute work of viscous stress and heat conduction
            Theta_ = (tau @ get_velocity(flow_face)) + flow_face["k"] * grad_T_

            # Compute viscous flux
            flux = zeros(5, dtype=float64)
            flux[1:4] = tau @ get_normal(face)
            flux[4] = Theta_.dot(get_normal(face))

            # Add the convective flux to the array of all convective fluxes
            Fv[left_cell_id, left_cell_face_idx, :] = -flux[:]
            Fv[right_cell_id, right_cell_face_idx, :] = flux[:]
