from numpy import array
from tqdm import tqdm
from time import time

from solver.Simulation import Simulation

def main():
    # Read mesh
    from solver.GmshGridReader import GmshGridReader
    mesh_reader = GmshGridReader()
    mesh_file = "/Users/maitreya/Desktop/pipe-3d-mesh/pipe.msh"
    nodes, faces, cells = mesh_reader.read_grid_from_file(mesh_file)

    # Initial conditions
    print("initializing flow field...")
    for cell in cells:
        # TEMP: Pipe Flow parabolic profile
        u_max =  1.
        R_pipe = 1.
        r = (cell.r_[1]**2 + cell.r_[2]**2)**0.5
        u = -u_max * (0.5 * (r / R_pipe)**2 + 0.5 * (r / R_pipe) - 1.)
        T = 273.
        p = 101325. # * (1. - 0.5 * (cell.r_[0] / 5.))
        rho = p / (287.058 * T)
        E = p / (0.4 * rho) + 0.5 * u**2
        cell.flow.W_ = rho * array([1., u, 0., 0., E])

    # Construct Simulation
    sim = Simulation(faces, cells)

    # TEMP: Time stepping
    num_iter = 10
    res_L2_norm = [None] * num_iter
    t = 0.
    start_time = time()
    for iter in range(num_iter):
        dt, res_L2_norm[iter] = sim.step(t)
        print("iteration {}, t = {}, dt = {}".format(iter, t, dt))
        t += dt
    end_time = time()
    iter_time = end_time - start_time
    print("total iteration time = {} sec".format(iter_time))
    print("time per iteration {} sec".format(iter_time / num_iter))

    sim.prepare_for_save()

    # TEMP: Debugging
    #from solver.SlipAdiabaticWallBC import SlipAdiabaticWallBC
    #from solver.InjectionBC import InjectionBC
    #from solver.OutletBC import OutletBC
    #import seaborn as sns
    #def cond(cell):
    #    return any([isinstance(face, BoundaryFace) and isinstance(face.bc, (OutletBC, InjectionBC)) for face in cell.faces]) and any([isinstance(face, BoundaryFace) and isinstance(face.bc, SlipAdiabaticWallBC) for face in cell.faces])
    #query = [dt*sum([cell.Fc_map[face] * face.area for face in cell.faces])[0]/cell.vol for cell in interior_cells if cond(cell)]
    #sns.distplot(np.array(query), rug=True, kde=False)
    import matplotlib.pyplot as plt
    plt.figure()
    for i in range(5):
        plt.semilogy(range(num_iter), [res[i] for res in res_L2_norm], label=i)
    plt.legend()
    plt.show()

    # Write data file
    from solver.VTKWriter import VTKWriter
    mesh_writer = VTKWriter()
    write_ghost = False
    output_file = "/Users/maitreya/Desktop/pipe-3d-mesh/pipe.vtk"
    mesh_writer.write_to_file(nodes, cells, output_file, write_ghost=write_ghost)

if __name__ == "__main__":
    main()
