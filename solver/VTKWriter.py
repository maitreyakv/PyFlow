import solver

# TODO: Add doc for class
class VTKWriter:

    cell_types = {solver.TetrahedronCell.TetrahedronCell: 10}

    # Constructor for VTKWriter
    def __init__(self):
        pass

    # TODO: Add doc for function
    def write_to_file(self, cells, filename):
        with open(filename, "w") as fp:
            fp.write("# vtk DataFile Version 2.0\ncomment goes here\nASCII\nDATASET UNSTRUCTURED_GRID\n\n")

            nodes = list({str(pt): pt for cell in cells for pt in cell.pts}.values())
            fp.write("POINTS {} double\n".format(len(nodes)))
            for i, node in enumerate(nodes):
                fp.write(" ".join(map(str, node)) + "\n")

            node_id_map = {str(node): id for id, node in enumerate(nodes)}
            def node_id(node):
                return node_id_map[str(node)]

            cell_str = ""
            cell_types_str = ""
            cell_section_size = 0
            for cell in cells:
                cell_str += "{} ".format(len(cell.pts)) + " ".join(map(str, list(map(node_id, cell.pts)))) + "\n"
                cell_section_size += len(cell.pts) + 1
                cell_types_str += "{}\n".format(self.cell_types[type(cell)])
            fp.write("\nCELLS {} {}\n".format(len(cells), cell_section_size))
            fp.write(cell_str)

            fp.write("\nCELL_TYPES {}\n".format(len(cells)))
            fp.write(cell_types_str)
