import numpy as np
from tqdm import tqdm

from solver.Simulation import Simulation
from solver.GhostCell import GhostCell
from solver.BoundaryFace import BoundaryFace

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
        T = 273.
        p = 101325. * (1. - 0.5 * (cell.r_[0] / 5.))
        rho = p / (287.058 * T)
        E = p / (0.4 * rho) + 0.5 * u**2
        cell.flow.W_ = rho * np.array([1., u, 0., 0., E])

    # Construct the simulation
    sim = Simulation(nodes, faces, cells)

    # Initialize Thermo
    from solver.CPGThermo import CPGThermo
    thermo = CPGThermo(1.4, 0.72, 287.058)

    # Apply BCs
    for face in faces:
        if isinstance(face, BoundaryFace):
            face.apply_bc_to_ghost_cell()

    # Update all variables in all cells
    for cell in cells:
        cell.flow.update(thermo)

    # Write data file
    from solver.VTKWriter import VTKWriter
    mesh_writer = VTKWriter()
    write_ghost = False
    mesh_writer.write_to_file(nodes, cells, "/Users/maitreya/Desktop/pipe-3d-mesh/pipe.vtk", write_ghost=write_ghost)

if __name__ == "__main__":
    main()
