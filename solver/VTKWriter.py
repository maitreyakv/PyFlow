from tqdm import tqdm

import solver

# TODO: Add doc for class
class VTKWriter:

    cell_types = {solver.TetrahedronCell.TetrahedronCell: 10}

    # Constructor for VTKWriter
    def __init__(self):
        pass

    # TODO: Add doc for function
    # TODO: Cleanup code and add comments
    def write_to_file(self, nodes, cells, filename):
        print("writing simulation data to VTK file...")

        with open(filename, "w") as fp:
            fp.write("# vtk DataFile Version 2.0\ncomment goes here\nASCII\nDATASET UNSTRUCTURED_GRID\n\n")

            print("writing mesh data...")

            fp.write("POINTS {} double\n".format(len(nodes)))
            for id in range(len(nodes)):
                fp.write(" ".join(map(str, nodes[id].r_.tolist())) + "\n")

            cell_str = ""
            cell_types_str = ""
            cell_section_size = 0
            for cell in tqdm(cells):
                cell_str += "{} ".format(len(cell.nodes)) + " ".join([str(node.id) for node in cell.nodes]) + "\n"
                cell_section_size += len(cell.nodes) + 1
                cell_types_str += "{}\n".format(self.cell_types[type(cell)])
            fp.write("\nCELLS {} {}\n".format(len(cells), cell_section_size))
            fp.write(cell_str)

            fp.write("\nCELL_TYPES {}\n".format(len(cells)))
            fp.write(cell_types_str)
