from numpy import zeros, float64, abs, sum, inf
from numpy.linalg import norm
from numba.core import types
from numba.typed import Dict

from solver.mesh import get_centroid, get_normal
from solver.flow import get_velocity
from solver.thermo import update_flow
from solver.bc import apply_bndry_cond
from solver.convective_flux import central_scheme_convective_fluxes
from solver.gradient import least_squares_gradient
from solver.viscous_flux import gradient_avg_viscous_fluxes

thermo = "cpg"
thermo_opts = Dict.empty(
    key_type=types.unicode_type,
    value_type=types.float64,
)
thermo_opts["gamma"] = 1.4
thermo_opts["Pr"] = 0.72
thermo_opts["R"] = 287.058

def compute_residual(faces, cells, flow_faces, flow_cells, W_faces, W_cells, t):
    # Update flow in alll cells
    update_flow(flow_cells, W_cells, thermo=thermo, opts=thermo_opts)

    # Set residuals to zero
    Rc = zeros((cells.size, 5), dtype=float64)
    Rd = zeros((cells.size, 5), dtype=float64)

    # Apply boundary conditions to boundary faces and ghost cells
    apply_bndry_cond(faces, cells, W_faces, W_cells, flow_cells, t, thermo=thermo, opts=thermo_opts)

    # Update flow in all cells
    update_flow(flow_cells, W_cells, thermo=thermo, opts=thermo_opts)

    # Allocate arrays for convective and viscous fluxes
    Fc = zeros((cells.size, 4, 5), dtype=float64)
    Fv = zeros((cells.size, 4, 5), dtype=float64)

    # Compute and store convective fluxes
    central_scheme_convective_fluxes(faces, cells, W_faces, W_cells, flow_faces, Fc, thermo="cpg", opts=thermo_opts)

    # TEMP: Debugging
    """
    #for cell in cells:
    #    if not cell["ghost"]:
    #        print( W_cells[cell["id"]][1:4] )
    #for face in faces:
    #    if face["bc"] >= 0:
    #        print(face["bc"], get_velocity(flow_faces[face["id"]]))

    import seaborn as sns
    import matplotlib.pyplot as plt
    from numpy import array, sum, any, all
    areas = [ array([ faces[cells[i]["face1"]]["area"],
                      faces[cells[i]["face2"]]["area"],
                      faces[cells[i]["face3"]]["area"],
                      faces[cells[i]["face4"]]["area"] ] ) for i in range(cells.size) ]
    areas = array(areas)
    mask1 = [ array([ faces[cells[i]["face1"]]["bc"] == 0,
                      faces[cells[i]["face2"]]["bc"] == 0,
                      faces[cells[i]["face3"]]["bc"] == 0,
                      faces[cells[i]["face4"]]["bc"] == 0 ] ) for i in range(cells.size) ]
    mask2 = [ array([ faces[cells[i]["face1"]]["bc"] not in [2,3],
                      faces[cells[i]["face2"]]["bc"] not in [2,3],
                      faces[cells[i]["face3"]]["bc"] not in [2,3],
                      faces[cells[i]["face4"]]["bc"] not in [2,3] ] ) for i in range(cells.size) ]
    #mask = any(array(mask1), axis=1) & ~all(array(mask2), axis=1) & ~cells["ghost"]
    mask = ~cells["ghost"]
    val = 7.7785e-05 * sum( Fc[:,:,:] * areas[:,:,None], axis=1 ) / cells["vol"][:,None]
    sns.distplot(val[mask,1], rug=True, kde=False)
    plt.show()
    exit()
    """

    # TODO: Add artificial dissipation

    # Compute gradients
    grads = least_squares_gradient(cells, flow_cells)

    # Compute and store viscous fluxes
    gradient_avg_viscous_fluxes(faces, cells, flow_faces, flow_cells, grads, Fv)

    # Add the convective and viscous fluxes to the residuals
    for cell in cells:
        if not cell["ghost"]:
            cell_id = cell["id"]
            area1 = faces[cell["face1"]]["area"]
            area2 = faces[cell["face2"]]["area"]
            area3 = faces[cell["face3"]]["area"]
            area4 = faces[cell["face4"]]["area"]
            Rc[cell_id,:] += Fc[cell_id,0,:] * area1 \
                           + Fc[cell_id,1,:] * area2 \
                           + Fc[cell_id,2,:] * area3 \
                           + Fc[cell_id,3,:] * area4
            Rd[cell_id,:] += Fv[cell_id,0,:] * area1 \
                           + Fv[cell_id,1,:] * area2 \
                           + Fv[cell_id,2,:] * area3 \
                           + Fv[cell_id,3,:] * area4

    # Return the residuals
    return Rc, Rd

# Five stage coefficients for Central Scheme
alpha1 = 0.25
alpha2 = 0.1667
alpha3 = 0.375
alpha4 = 0.5
alpha5 = 1.
beta3 = 0.56
beta5 = 0.44

# CFL number for central scheme
CFL = 3.6

# Spectral radii constant for central scheme
C = 4.

def hybrid_multi_stage_integrate(faces, cells, flow_faces, flow_cells, W_faces, W_cells, t):
    # Compute area projections when the function runs for the first time
    if type(hybrid_multi_stage_integrate.area_proj) == type(None):
        area_proj = zeros((cells.size, 3), dtype=float64)
        for cell in cells:
            if not cell["ghost"]:
                face1 = faces[cell["face1"]]
                face2 = faces[cell["face2"]]
                face3 = faces[cell["face3"]]
                face4 = faces[cell["face4"]]
                r_cell_ = get_centroid(cell)
                area1_ = get_normal(face1, ref_point=r_cell_) * face1["area"]
                area2_ = get_normal(face2, ref_point=r_cell_) * face2["area"]
                area3_ = get_normal(face3, ref_point=r_cell_) * face3["area"]
                area4_ = get_normal(face4, ref_point=r_cell_) * face4["area"]
                area_proj[cell["id"],:] = 0.5 * ( abs(area1_) + abs(area2_) + abs(area3_) + abs(area4_) )
        hybrid_multi_stage_integrate.area_proj = area_proj

    # Get stored area projections for the cells
    area_proj = hybrid_multi_stage_integrate.area_proj

    # Compute initial residual
    Rc0, Rd0 = compute_residual(faces, cells, flow_faces, flow_cells, W_faces, W_cells, t)

    # Compute time step
    dt = inf
    for cell in cells:
        if not cell["ghost"]:
            cell_id = cell["id"]
            flow_cell = flow_cells[cell_id]
            spec_rad_c = (abs(get_velocity(flow_cell)) + flow_cell["c"]) * area_proj[cell_id,:]
            spec_rad_v = max(4. / (3. * flow_cell["rho"]), flow_cell["gamma"] / flow_cell["rho"]) * (flow_cell["k"] / flow_cell["cp"]) * area_proj[cell_id]**2 / cell["vol"]
            dt_local = CFL * cell["vol"] / (sum(spec_rad_c) + C * sum(spec_rad_v))
            dt = min(dt, dt_local)

    # Perform first stage
    W1 = W_cells - alpha1 * dt * (Rc0 - Rd0) / cells["vol"][:,None]

    # Perform second stage
    Rc1, _ = compute_residual(faces, cells, flow_faces, flow_cells, W_faces, W1, t)
    W2 = W_cells - alpha2 * dt * (Rc1 - Rd0) / cells["vol"][:,None]

    # Perform third stage
    Rc2, Rd2 = compute_residual(faces, cells, flow_faces, flow_cells, W_faces, W2, t)
    Rd20 = beta3 * Rd2 + (1. - beta3)*Rd0
    W3 = W_cells - alpha3 * dt * (Rc2 - Rd20) / cells["vol"][:,None]

    # Perform fourth stage
    Rc3, _ = compute_residual(faces, cells, flow_faces, flow_cells, W_faces, W3, t)
    W4 = W_cells - alpha4 * dt * (Rc3 - Rd20) / cells["vol"][:,None]

    # Perform fifth stage
    Rc4, Rd4 = compute_residual(faces, cells, flow_faces, flow_cells, W_faces, W4, t)
    Rd42 = beta5 * Rd4 + (1. - beta5)*Rd20
    W_cells[:,:] = W_cells - alpha5 * dt * (Rc4 - Rd42) / cells["vol"][:,None]

    # Compute residual of last stage and L2 norm
    # TEMP: Debugging
    res = dt * (Rc0 - Rd0) / cells["vol"][:,None] #Rc4 - Rd42
    res_L2_norm = norm(res, axis=0)

    # Perform final update of flow in all cells
    update_flow(flow_cells, W_cells, thermo=thermo, opts=thermo_opts)

    return dt, res_L2_norm

hybrid_multi_stage_integrate.area_proj = None
