from solver.mesh import get_normal, get_centroid
from solver.thermo import E, rho

def apply_bndry_cond(faces, cells, W_faces, W_cells, flow_cells, t, thermo="cpg", opts=None):
    for face in faces:
        if face["bc"] > 0:
            inter_cell = cells[face["left_cell"]]
            ghost_cell = cells[face["left_cell"]]
            W_bndry_face = W_faces[face["id"]]
            W_inter_cell = W_cells[inter_cell["id"]]
            W_ghost_cell = W_cells[ghost_cell["id"]]
            flow_inter_cell = flow_cells[inter_cell["id"]]

            if face["bc"] == 1: # TEMP: Wall bc implementation
                slip_wall_bc(face, inter_cell, ghost_cell, W_bndry_face, W_inter_cell,
                             W_ghost_cell, flow_inter_cell, t, thermo=thermo, opts=opts)
            elif face["bc"] == 2: # TEMP: Outlet bc implementation
                outlet_bc(face, inter_cell, ghost_cell, W_bndry_face, W_inter_cell,
                          W_ghost_cell, flow_inter_cell, t, thermo=thermo, opts=opts)
            elif face["bc"] == 3: # TEMP: Injection bc implementation
                injection_bc(face, inter_cell, ghost_cell, W_bndry_face, W_inter_cell,
                             W_ghost_cell, flow_inter_cell, t, thermo=thermo, opts=opts)
            else:
                print("error: cannot apply boundary condition")

def slip_wall_bc(bndry_face, inter_cell, ghost_cell, W_bndry_face, W_inter_cell, W_ghost_cell,
                    flow_inter_cell, t, thermo="cpg", opts=None):
    # Copy the density and energy from the interior cell to the ghost cell
    W_ghost_cell[0] = W_inter_cell[0]
    W_ghost_cell[4] = W_inter_cell[4]

    # Compute interior cell velocity
    v_inter_ = W_inter_cell[1:4] / W_inter_cell[0]

    # Compute boundary face normal vector
    n_ = get_normal(bndry_face)

    # Compute contravariant velocity
    V = v_inter_.dot(n_)

    # Reflect velocity from the interior cell to the ghost cell
    W_ghost_cell[1:4] = W_ghost_cell[0] * (v_inter_ - 2. * V * n_)

    # Interpolate conservative variables to the BoundaryFace
    W_bndry_face[:] = 0.5 * (W_ghost_cell[:] + W_inter_cell[:])

def outlet_bc(bndry_face, inter_cell, ghost_cell, W_bndry_face, W_inter_cell, W_ghost_cell,
                flow_inter_cell, t, thermo="cpg", opts=None):
    # Define the reference state using the interior cell
    rho_ref = flow_inter_cell["rho"]
    c_ref = flow_inter_cell["c"]

    # Use the outflow pressure as the boundary pressure
    # TEMP: Hardcoded pressure
    p_bndry = 101325.

    # Get interior cell flow velocity
    v_inter_ = W_inter_cell[1:4] / W_inter_cell[0]

    # Compute boundary face normal vector
    n_ = get_normal(bndry_face)

    # Compute boundary flow properties
    rho_bndry = flow_inter_cell["rho"] + (p_bndry - flow_inter_cell["p"]) / c_ref**2
    v_bndry_ = v_inter_ + n_ * (flow_inter_cell["p"] - p_bndry) / (rho_ref * c_ref)

    # Compute the energy on the boundary
    E_bndry = E(p_bndry, rho_bndry, v_bndry_, thermo=thermo, opts=opts)

    # Update conservative variables in the boundary
    W_bndry_face[0] = rho_bndry
    W_bndry_face[1:4] = rho_bndry * v_bndry_
    W_bndry_face[4] = rho_bndry * E_bndry

    # Extrapolate conservative variables in the ghost cell from the boundary and interior
    W_ghost_cell[:] = 2. * W_bndry_face[:] - W_inter_cell[:]

def injection_bc(bndry_face, inter_cell, ghost_cell, W_bndry_face, W_inter_cell, W_ghost_cell,
                    flow_inter_cell, t, thermo="cpg", opts=None):
    # Set boundary pressure as interior cell pressure
    p_bndry = flow_inter_cell["p"]

    # Set boundary temperature as injection temperature
    # TEMP: Hardocded temperature
    T_bndry = 273.

    # Compute boundary density
    rho_bndry = rho(p_bndry, T_bndry, thermo=thermo, opts=opts)

    # Compute boundary face normal vector
    n_ = get_normal(bndry_face)

    # Compute boundary velocity
    # TEMP: Hardcoded mass flow/velocity
    v_bndry_ = -n_ * 1.

    # Compute boundary energy
    E_bndry = E(p_bndry, rho_bndry, v_bndry_, thermo=thermo, opts=opts)

    # Set the boundary conservative flow variables
    W_bndry_face[0] = rho_bndry
    W_bndry_face[1:4] = rho_bndry * v_bndry_
    W_bndry_face[4] = rho_bndry * E_bndry

    # Extrapolate conservative variables to ghost cell
    W_ghost_cell[:] = 2. * W_bndry_face[:] - W_inter_cell[:]
