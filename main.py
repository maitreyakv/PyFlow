import numpy as np
from tqdm import tqdm

from solver.Simulation import Simulation

def main():
    # Read mesh
    from solver.GmshGridReader import GmshGridReader
    mesh_reader = GmshGridReader()
    nodes, faces, cells = mesh_reader.read_grid_from_file("/Users/maitreya/Desktop/pipe-3d-mesh/pipe.msh")

    # Initial conditions
    print("initializing flow field...")
    for cell in cells:
        # TEMP: Pipe Flow parabolic profile
        u_max = 1.
        R_pipe = 1.
        r = np.sqrt(cell.r_[1]**2 + cell.r_[2]**2)
        u = -u_max * (0.5 * (r / R_pipe)**2 + 0.5 * (r / R_pipe) - 1.)
        cell.flow.W_ = 1. * np.array([1., u, 0., 0., 1.])

    # Construct the simulation
    sim = Simulation(nodes, faces, cells)

    # Write data file
    from solver.VTKWriter import VTKWriter
    mesh_writer = VTKWriter()
    mesh_writer.write_to_file(nodes, cells, "/Users/maitreya/Desktop/pipe-3d-mesh/pipe.vtk")

if __name__ == "__main__":
    main()
