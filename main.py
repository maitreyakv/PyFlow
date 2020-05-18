from numpy import array
from tqdm import tqdm

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
        u = u_max # -u_max * (0.5 * (r / R_pipe)**2 + 0.5 * (r / R_pipe) - 1.)
        T = 273.
        p = 101325. # * (1. - 0.5 * (cell.r_[0] / 5.))
        rho = p / (287.058 * T)
        E = p / (0.4 * rho) + 0.5 * u**2
        cell.flow.W_ = rho * array([1., u, 0., 0., E])

    # Construct Simulation
    sim = Simulation(faces, cells)

    # TEMP: Time stepping
    t = 0.
    dt = 0.01 * 0.3 / 1.
    for iter in range(1):
        print("iteration {}...".format(iter))

        sim.compute_residuals(t)

        sim.integrate_step(dt)

        t += dt

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
    #import matplotlib.pyplot as plt
    #plt.show()

    # Write data file
    from solver.VTKWriter import VTKWriter
    mesh_writer = VTKWriter()
    write_ghost = False
    output_file = "/Users/maitreya/Desktop/pipe-3d-mesh/pipe.vtk"
    mesh_writer.write_to_file(nodes, cells, output_file, write_ghost=write_ghost)

if __name__ == "__main__":
    main()
