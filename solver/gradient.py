from numpy import zeros, float64, array
from numpy.linalg import lstsq, qr, inv

from solver.mesh import get_centroid

# TODO: Add doc and cleanup

def least_squares_gradient(cells, flow_cells):
    # Compute edge matrices for each cell when function is called for the first time
    if type(least_squares_gradient.edge_matrices) == type(None):
        edge_matrices = zeros((cells.size, 4, 3), dtype=float64)
        R_inv_Q_t = zeros((cells.size, 3, 4), dtype=float64)
        for cell in cells:
            if not cell["ghost"]:
                r_cell_ = get_centroid(cell)
                edge_matrices[cell["id"],0,:] = get_centroid(cells[cell["nbr1"]]) - r_cell_
                edge_matrices[cell["id"],1,:] = get_centroid(cells[cell["nbr2"]]) - r_cell_
                edge_matrices[cell["id"],2,:] = get_centroid(cells[cell["nbr3"]]) - r_cell_
                edge_matrices[cell["id"],3,:] = get_centroid(cells[cell["nbr4"]]) - r_cell_
                Q, R = qr(edge_matrices[cell["id"],:,:])
                R_inv_Q_t[cell["id"],:,:] = inv(R) @ Q.T
        least_squares_gradient.edge_matrices = edge_matrices
        least_squares_gradient.R_inv_Q_t = R_inv_Q_t

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

            #edge_matrix = least_squares_gradient.edge_matrices[cell["id"]]
            R_inv_Q_t = least_squares_gradient.R_inv_Q_t[cell["id"],:,:]

            grads[cell["id"],0,:] = R_inv_Q_t @ diff_u #lstsq(edge_matrix, diff_u, rcond=None)[0]
            grads[cell["id"],1,:] = R_inv_Q_t @ diff_v #lstsq(edge_matrix, diff_v, rcond=None)[0]
            grads[cell["id"],2,:] = R_inv_Q_t @ diff_w #lstsq(edge_matrix, diff_w, rcond=None)[0]
            grads[cell["id"],3,:] = R_inv_Q_t @ diff_T #lstsq(edge_matrix, diff_T, rcond=None)[0]

    # Return the computed gradients
    return grads

least_squares_gradient.edge_matrices = None
