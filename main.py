from numpy import zeros, float64, array
from time import time

from solver.io import read_gmsh_2, write_vtk
from solver.flow import flow_type
from solver.integration import hybrid_multi_stage_integrate

def main():
    nodes, faces, cells = read_gmsh_2("/Users/maitreya/Desktop/pipe-3d-mesh/pipe.msh")

    flow_cells = zeros(cells.size, dtype=flow_type)
    flow_faces = zeros(faces.size, dtype=flow_type)

    num_cons_var = 5
    W_cells = zeros((cells.size, num_cons_var), dtype=float64)
    W_faces = zeros((faces.size, num_cons_var), dtype=float64)

    # TEMP: Pipe flow initial conditions
    u = 1.
    T = 273.
    p = 101325.
    rho = p / (287.058 * T)
    E = p / (0.4 * rho) + 0.5 * u**2
    W_init = rho * array( [ 1., u, 0., 0., E ] )
    W_cells[:,:] = W_init[None,:]
    W_faces[:,:] = W_init[None,:]

    t = 0.
    num_iter = 4
    start_time = time()
    for iter in range(num_iter):
        dt = hybrid_multi_stage_integrate(faces, cells, flow_faces, flow_cells, W_faces, W_cells, t)
        t += dt
        print("iter {}, t={}, dt={}".format(iter, t, dt))
    end_time = time()
    iter_time = end_time - start_time
    print("total iteration time = {} sec".format(iter_time))
    print("time per iteration {} sec".format(iter_time / num_iter))

    write_vtk("/Users/maitreya/Desktop/pipe-3d-mesh/pipe.vtk", nodes, faces, cells, ghost=False)

if __name__ == "__main__":
    main()
