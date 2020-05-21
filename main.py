from numpy import zeros, float64, array, array2string
from time import time

from solver.io import read_gmsh_2, write_vtk
from solver.flow import flow_type
from solver.integration import hybrid_multi_stage_integrate

def main():
    gmsh_filename = "/Users/maitreya/Desktop/channel-3d-mesh/channel.msh"
    nodes, faces, cells = read_gmsh_2(gmsh_filename)

    flow_cells = zeros(cells.size, dtype=flow_type)
    flow_faces = zeros(faces.size, dtype=flow_type)

    num_cons_var = 5
    W_cells = zeros((cells.size, num_cons_var), dtype=float64)
    W_faces = zeros((faces.size, num_cons_var), dtype=float64)

    # TEMP: Pipe flow initial conditions
    u = 100.
    T = 273.
    p = 101325.
    rho = p / (287.058 * T)
    E = p / (0.4 * rho) + 0.5 * u**2
    W_init = rho * array( [ 1., u, 0., 0., E ] )
    W_cells[:,:] = W_init[None,:]
    W_faces[:,:] = W_init[None,:]

    t = 0.
    num_iter = 20
    res_L2_norm_hist = zeros((num_iter, 5), dtype=float64)
    start_time = time()
    for iter in range(num_iter):
        dt, res_L2_norm = hybrid_multi_stage_integrate(faces, cells, flow_faces, flow_cells,
                                                                            W_faces, W_cells, t)
        t += dt
        res_L2_norm_hist[iter,:] = res_L2_norm[:]
        print("iter {:4d}, t={:.4e}, dt={:.4e}, L2 norm of res.={}".format(
              iter, t, dt, array2string(res_L2_norm).replace('\n', '')))
    end_time = time()
    iter_time = end_time - start_time
    print("total iteration time = {} sec".format(iter_time))
    print("time per iteration {} sec".format(iter_time / num_iter))

    # Plot history of residuals
    from solver.plot import plot_residual_history
    plot_residual_history(res_L2_norm_hist)

    # Write final data to output file
    vtk_filename = "/Users/maitreya/Desktop/channel-3d-mesh/channel.vtk"
    write_vtk(vtk_filename, nodes, faces, cells, flow_cells, ghost=False)

if __name__ == "__main__":
    main()
