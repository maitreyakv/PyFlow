from numpy import zeros, float64
from numba.core import types
from numba.typed import Dict

from solver.io import read_gmsh_2, write_vtk
from solver.flow import flow_type
from solver.thermo import update_flow

def main():
    nodes, faces, cells = read_gmsh_2("/Users/maitreya/Desktop/pipe-3d-mesh/pipe.msh")

    flow_cells = zeros(cells.size, dtype=flow_type)
    flow_faces = zeros(faces.size, dtype=flow_type)

    num_cons_var = 5
    W_cells = zeros((cells.size, num_cons_var), dtype=float64)
    W_faces = zeros((cells.size, num_cons_var), dtype=float64)

    thermo = "cpg"
    thermo_opts = Dict.empty(
        key_type=types.unicode_type,
        value_type=types.float64,
    )
    thermo_opts["gamma"] = 1.4
    thermo_opts["Pr"] = 0.72
    thermo_opts["R"] = 287.058

    update_flow(flow_cells, W_cells, thermo=thermo, opts=thermo_opts)

    write_vtk("/Users/maitreya/Desktop/pipe-3d-mesh/pipe.vtk", nodes, faces, cells)

if __name__ == "__main__":
    main()
