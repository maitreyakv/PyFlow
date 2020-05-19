from numpy import zeros, float64, array, ascontiguousarray, asfortranarray
from numpy.linalg import lstsq, qr, inv
from numba import jit

from solver.mesh import get_centroid

# TODO: Add doc and cleanup

@jit(nopython=True)
def least_squares_gradient(cells, flow_cells):
    # Compute edge matrices for each cell
    QR_fac = zeros((cells.size, 3, 4), dtype=float64)
    for cell in cells:
        if not cell["ghost"]:
            edge_matrix = zeros((4, 3), dtype=float64)
            r_cell_ = get_centroid(cell)
            edge_matrix[0,:] = get_centroid(cells[cell["nbr1"]]) - r_cell_
            edge_matrix[1,:] = get_centroid(cells[cell["nbr2"]]) - r_cell_
            edge_matrix[2,:] = get_centroid(cells[cell["nbr3"]]) - r_cell_
            edge_matrix[3,:] = get_centroid(cells[cell["nbr4"]]) - r_cell_
            Q, R = qr(edge_matrix)
            QR_fac[cell["id"],:,:] = ascontiguousarray( inv(R) ) @ ascontiguousarray( Q.T )

    # Allocate arrays for the gradients
    grads = zeros((cells.size, 4, 3), dtype=float64)

    # Compute gradients in each cell using least square gradients
    for cell in cells:
        if not cell["ghost"]:
            flow_cell = flow_cells[cell["id"]]
            flow_nbr1 = flow_cells[cell["nbr1"]]
            flow_nbr2 = flow_cells[cell["nbr2"]]
            flow_nbr3 = flow_cells[cell["nbr3"]]
            flow_nbr4 = flow_cells[cell["nbr4"]]
            diff_u = array( [ flow_nbr1["u"] - flow_cell["u"],
                              flow_nbr2["u"] - flow_cell["u"],
                              flow_nbr3["u"] - flow_cell["u"],
                              flow_nbr4["u"] - flow_cell["u"] ] )
            diff_v = array( [ flow_nbr1["v"] - flow_cell["v"],
                              flow_nbr2["v"] - flow_cell["v"],
                              flow_nbr3["v"] - flow_cell["v"],
                              flow_nbr4["v"] - flow_cell["v"] ] )
            diff_w = array( [ flow_nbr1["w"] - flow_cell["w"],
                              flow_nbr2["w"] - flow_cell["w"],
                              flow_nbr3["w"] - flow_cell["w"],
                              flow_nbr4["w"] - flow_cell["w"] ] )
            diff_T = array( [ flow_nbr1["T"] - flow_cell["T"],
                              flow_nbr2["T"] - flow_cell["T"],
                              flow_nbr3["T"] - flow_cell["T"],
                              flow_nbr4["T"] - flow_cell["T"] ] )

            #R_inv_Q_t = zeros((3, 4), dtype=float64)
            R_inv_Q_t = ascontiguousarray( QR_fac[cell["id"],:,:] )

            grads[cell["id"],0,:] = R_inv_Q_t @ diff_u
            grads[cell["id"],1,:] = R_inv_Q_t @ diff_v
            grads[cell["id"],2,:] = R_inv_Q_t @ diff_w
            grads[cell["id"],3,:] = R_inv_Q_t @ diff_T

    # Return the computed gradients
    return grads
