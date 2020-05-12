from solver.Cell import Cell
from solver.GhostCell import GhostCell
from solver.TetrahedronCell import TetrahedronCell

# TODO: Add doc for class
class VTKWriter:

    cell_types = {TetrahedronCell: 10,
                  GhostCell: 10}

    # Constructor for VTKWriter
    def __init__(self):
        pass

    # TODO: Add doc for function
    # TODO: Cleanup code and add comments
    def write_to_file(self, nodes, cells, filename, write_ghost=False):
        print("writing simulation data to VTK file...")

        skip_cell_type = NoneType if write_ghost else GhostCell


        with open(filename, "w") as fp:
            fp.write("# vtk DataFile Version 2.0\ncomment goes here\nASCII\nDATASET UNSTRUCTURED_GRID\n\n")

            print("writing mesh data...")

            fp.write("POINTS {} double\n".format(len(nodes)))
            for id in range(len(nodes)):
                fp.write(" ".join(map(str, nodes[id].r_.tolist())) + "\n")

            num_cells = len([cell for cell in cells if not isinstance(cell, skip_cell_type)])

            cell_str = ""
            cell_types_str = ""
            cell_section_size = 0
            for cell in cells:
                if not isinstance(cell, skip_cell_type):
                    cell_str += "{} ".format(len(cell.nodes)) + " ".join([str(node.id) for node in cell.nodes]) + "\n"
                    cell_section_size += len(cell.nodes) + 1
                    cell_types_str += "{}\n".format(self.cell_types[type(cell)])
            fp.write("\nCELLS {} {}\n".format(num_cells, cell_section_size))
            fp.write(cell_str)

            fp.write("\nCELL_TYPES {}\n".format(num_cells))
            fp.write(cell_types_str)

            print("writing flow data...")
            fp.write("CELL_DATA {}\n".format(num_cells))

            fp.write("SCALARS pressure double 1\nLOOKUP_TABLE default\n")
            for cell in cells:
                if not isinstance(cell, skip_cell_type):
                    fp.write("{}\n".format(cell.flow.p))

            fp.write("\nSCALARS density double 1\nLOOKUP_TABLE default\n")
            for cell in cells:
                if not isinstance(cell, skip_cell_type):
                    fp.write("{}\n".format(cell.flow.rho))

            fp.write("\nSCALARS temperature double 1\nLOOKUP_TABLE default\n")
            for cell in cells:
                if not isinstance(cell, skip_cell_type):
                    fp.write("{}\n".format(cell.flow.T))

            fp.write("\nVECTORS velocity double\n")
            for cell in cells:
                if not isinstance(cell, skip_cell_type):
                    fp.write(" ".join(map(str, cell.flow.v_)) + '\n')
