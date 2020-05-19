from numpy import zeros, float64
from numba.core import types
from numba.typed import Dict

from solver.thermo import update_flow
from solver.bc import apply_bndry_cond

thermo = "cpg"
thermo_opts = Dict.empty(
    key_type=types.unicode_type,
    value_type=types.float64,
)
thermo_opts["gamma"] = 1.4
thermo_opts["Pr"] = 0.72
thermo_opts["R"] = 287.058

def compute_residual(faces, cells, flow_faces, flow_cells, W_faces, W_cells, t):
    # Update flow in cells
    update_flow(flow_cells, W_cells, thermo=thermo, opts=thermo_opts)

    # Set residuals to zero
    Rc = zeros((cells.size, 5), dtype=float64)
    Rd = zeros((cells.size, 5), dtype=float64)

    # Apply boundary conditions to boundary faces and ghost cells
    apply_bndry_cond(faces, cells, W_faces, W_cells, flow_cells, t, thermo=thermo, opts=thermo_opts)

    # Return the residuals
    return Rc, Rd


def hybrid_multi_stage_integrate(faces, cells, flow_faces, flow_cells, W_faces, W_cells, t):
    Rc, Rd = compute_residual(faces, cells, flow_faces, flow_cells, W_faces, W_cells, t)
