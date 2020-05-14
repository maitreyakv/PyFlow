import numpy as np
from tqdm import tqdm

from solver.GhostCell import GhostCell
from solver.BoundaryFace import BoundaryFace

def main():
    # Read mesh
    from solver.GmshGridReader import GmshGridReader
    mesh_reader = GmshGridReader()
    nodes, faces, cells = mesh_reader.read_grid_from_file("/Users/maitreya/Desktop/pipe-3d-mesh/pipe.msh")

    interior_cells = [cell for cell in cells if not isinstance(cell, GhostCell)]
    ghost_cells = [cell for cell in cells if isinstance(cell, GhostCell)]

    boundary_faces = [face for face in faces if isinstance(face, BoundaryFace)]

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

    # Initialize Thermo
    from solver.CPGThermo import CPGThermo
    thermo = CPGThermo(1.4, 0.72, 287.058)

    # Update all varaibles in interior cells
    for interior_cell in interior_cells:
        interior_cell.flow.update(thermo)

    # Apply BCs to GhostCells and BoundaryFaces
    for boundary_face in boundary_faces:
        boundary_face.apply_bc(thermo)

    # Update all variables in GhostCells
    for ghost_cell in ghost_cells:
        ghost_cell.flow.update(thermo)

    # Write data file
    from solver.VTKWriter import VTKWriter
    mesh_writer = VTKWriter()
    write_ghost = True
    mesh_writer.write_to_file(nodes, cells, "/Users/maitreya/Desktop/pipe-3d-mesh/pipe.vtk", write_ghost=write_ghost)

if __name__ == "__main__":
    main()
