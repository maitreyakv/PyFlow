import numpy as np
from tqdm import tqdm

from solver.GhostCell import GhostCell
from solver.BoundaryFace import BoundaryFace

def main():
    # Read mesh
    from solver.GmshGridReader import GmshGridReader
    mesh_reader = GmshGridReader()
    mesh_file = "/Users/maitreya/Desktop/pipe-3d-mesh/pipe.msh"
    nodes, faces, cells = mesh_reader.read_grid_from_file(mesh_file)

    interior_cells = [cell for cell in cells if not isinstance(cell, GhostCell)]
    ghost_cells = [cell for cell in cells if isinstance(cell, GhostCell)]

    boundary_faces = [face for face in faces if isinstance(face, BoundaryFace)]

    # Initial conditions
    print("initializing flow field...")
    for cell in cells:
        # TEMP: Pipe Flow parabolic profile
        u_max = 0. # 1.
        R_pipe = 1.
        r = np.sqrt(cell.r_[1]**2 + cell.r_[2]**2)
        u = -u_max * (0.5 * (r / R_pipe)**2 + 0.5 * (r / R_pipe) - 1.)
        T = 273.
        p = 101325. # * (1. - 0.5 * (cell.r_[0] / 5.))
        rho = p / (287.058 * T)
        E = p / (0.4 * rho) + 0.5 * u**2
        cell.flow.W_ = rho * np.array([1., u, 0., 0., E])

    # Initialize Thermo
    from solver.CPGThermo import CPGThermo
    thermo = CPGThermo(1.4, 0.72, 287.058)

    # Initialize Convective Flux Scheme
    from solver.CentralConvectiveFlux import CentralConvectiveFlux
    convective_flux = CentralConvectiveFlux(interior_cells)

    # Initialise Gradient computer
    from solver.GreenGaussGradient import GreenGaussGradient
    gradient = GreenGaussGradient()

    # Initialize Viscous Flux Scheme
    from solver.GradientAvgViscousFlux import GradientAvgViscousFlux
    viscous_flux = GradientAvgViscousFlux()

    # TEMP: Time stepping
    t = 0.
    dt = 0.1 * 0.25 / 1.
    for iter in range(2):
        print("iteration {}...".format(iter))

        # Update all variables in interior cells and set residuals to zero
        for interior_cell in interior_cells:
            interior_cell.flow.update(thermo)
            interior_cell.residual = np.zeros(5)

        # Apply BCs to GhostCells and BoundaryFaces and update all variables on BoundaryFaces
        for boundary_face in boundary_faces:
            boundary_face.apply_bc(thermo, t)
            boundary_face.flow.update(thermo)

        # Update all variables in GhostCells
        for ghost_cell in ghost_cells:
            ghost_cell.flow.update(thermo)

        # Compute the convective fluxes
        #print("computing convective fluxes...")
        for face in faces:
            convective_flux.compute_convective_flux(face, thermo)

        # Update Cell's spectral radius of convective flux jacobian
        for cell in interior_cells:
            convective_flux.prepare_for_artificial_dissipation(cell)

        # Add the artificial dissipation
        for cell in interior_cells:
            convective_flux.add_artificial_dissipation(cell)

        # Compute velocity and temperature gradients
        for cell in interior_cells:
            gradient.compute_gradients(cell)

        # Compute the viscous fluxes
        #print("computing viscous fluxes...")
        for face in faces:
            viscous_flux.compute_viscous_flux(face, thermo)

        # Add the fluxes to the residuals and compute L2 norm of residuals
        residual_L2_norm = np.zeros(5)
        for cell in interior_cells:
            cell.add_fluxes_to_residual()
            residual_L2_norm += cell.residual**2
        residual_L2_norm = np.sqrt(residual_L2_norm)
        print("residual L2 norm = {}".format(residual_L2_norm))

        # TEMP: Perform Euler integration
        for cell in interior_cells:
            cell.flow.W_ += -dt * cell.residual / cell.vol
        t += dt

    # Update all variables in all cells before writing to file
    for cell in cells:
        cell.flow.update(thermo)

    # TEMP: Debugging
    from solver.NoSlipAdiabaticWallBC import NoSlipAdiabaticWallBC
    #from solver.InjectionBC import InjectionBC
    #from solver.OutletBC import OutletBC
    import seaborn as sns
    def cond(cell):
        return True
    #    #return any([isinstance(face, BoundaryFace) and isinstance(face.bc, NoSlipAdiabaticWallBC) for face in cell.faces])
    query_cells = [sum(cell.Fc_map.values())[1] for cell in interior_cells if cond(cell)]
    sns.distplot(np.array(query_cells), rug=True)
    import matplotlib.pyplot as plt
    plt.show()

    # Write data file
    from solver.VTKWriter import VTKWriter
    mesh_writer = VTKWriter()
    write_ghost = False
    output_file = "/Users/maitreya/Desktop/pipe-3d-mesh/pipe.vtk"
    mesh_writer.write_to_file(nodes, cells, output_file, write_ghost=write_ghost)

if __name__ == "__main__":
    main()
