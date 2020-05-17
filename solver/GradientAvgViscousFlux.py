import numpy as np

from solver.ViscousFlux import ViscousFlux
from solver.BoundaryFace import BoundaryFace

# TODO: Add doc for class
class GradientAvgViscousFlux(ViscousFlux):

    # Constructor for GradientAvgViscousFlux
    def __init__(self):
        # Superclass constructor
        super().__init__()

    # Implements the computation of the viscous flux at a Face
    def compute_viscous_flux(self, face, thermo):
        # If BoundaryFace, then use interior Cell gradients as BoundaryFace gradients
        if isinstance(face, BoundaryFace):
            face.flow.grad_u_ = face.left_cell.flow.grad_u_
            face.flow.grad_v_ = face.left_cell.flow.grad_v_
            face.flow.grad_w_ = face.left_cell.flow.grad_w_
            face.flow.grad_T_ = face.left_cell.flow.grad_T_
        else:
            # Compute simple average of gradients at left and right cell
            grad_u_avg_ = 0.5 * (face.left_cell.flow.grad_u_ + face.right_cell.flow.grad_u_)
            grad_v_avg_ = 0.5 * (face.left_cell.flow.grad_v_ + face.right_cell.flow.grad_v_)
            grad_w_avg_ = 0.5 * (face.left_cell.flow.grad_w_ + face.right_cell.flow.grad_w_)
            grad_T_avg_ = 0.5 * (face.left_cell.flow.grad_T_ + face.right_cell.flow.grad_T_)

            # Compute distance between centroids of left and right Cells
            l = np.linalg.norm(face.right_cell.r_ - face.left_cell.r_)

            # Compute unit vector along this line
            t_ = (face.right_cell.r_ - face.left_cell.r_) / l

            # Compute directional derivative along line connecting left and right Cell centroids
            dudl = (face.right_cell.flow.v_[0] - face.left_cell.flow.v_[0]) / l
            dvdl = (face.right_cell.flow.v_[1] - face.left_cell.flow.v_[1]) / l
            dwdl = (face.right_cell.flow.v_[2] - face.left_cell.flow.v_[2]) / l
            dTdl = (face.right_cell.flow.T - face.left_cell.flow.T) / l

            # Compute modified average of gradients
            face.flow.grad_u_ = grad_u_avg_ - (grad_u_avg_.dot(t_) - dudl) * t_
            face.flow.grad_v_ = grad_v_avg_ - (grad_v_avg_.dot(t_) - dvdl) * t_
            face.flow.grad_w_ = grad_w_avg_ - (grad_w_avg_.dot(t_) - dwdl) * t_
            face.flow.grad_T_ = grad_T_avg_ - (grad_T_avg_.dot(t_) - dTdl) * t_

        # Compute divergence of velocity
        div_ = face.flow.grad_u_[0] + face.flow.grad_v_[1] + face.flow.grad_w_[2]

        # Assemble stress tensor
        tau = np.zeros((3,3))
        tau[0,0] = 2. * face.flow.mu * (face.flow.grad_u_[0] - div_ / 3.)
        tau[1,1] = 2. * face.flow.mu * (face.flow.grad_v_[1] - div_ / 3.)
        tau[2,2] = 2. * face.flow.mu * (face.flow.grad_w_[2] - div_ / 3.)
        tau[0,1] = tau[1,0] = face.flow.mu * (face.flow.grad_u_[1] + face.flow.grad_v_[0])
        tau[0,2] = tau[2,0] = face.flow.mu * (face.flow.grad_u_[2] + face.flow.grad_w_[0])
        tau[1,2] = tau[2,1] = face.flow.mu * (face.flow.grad_v_[2] + face.flow.grad_w_[1])

        # Compute work of viscous stress and heat conduction
        Theta_ = (tau @ face.flow.v_) + face.flow.k * face.flow.grad_T_

        # Compute viscous flux
        Fv = np.zeros(5)
        Fv[1:4] = tau @ face.n_
        Fv[4] = Theta_.dot(face.n_)

        # Add the viscous flux to the map of viscous fluxes in the left and right Cell
        face.left_cell.Fv_map[face] = -Fv
        face.right_cell.Fv_map[face] = Fv
