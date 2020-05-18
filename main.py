from solver.io import read_gmsh_2, write_vtk

def main():
    nodes, faces, cells = read_gmsh_2("/Users/maitreya/Desktop/pipe-3d-mesh/pipe.msh")

    write_vtk("/Users/maitreya/Desktop/pipe-3d-mesh/pipe.vtk", nodes, faces, cells)

if __name__ == "__main__":
    main()
