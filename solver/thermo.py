from numpy import sqrt
from numpy.linalg import norm

from numba import jit

# TODO: Add doc and cleanup

@jit(nopython=True)
def update_flow(flow, W, thermo="cpg", opts=None):
    flow["rho"][:] = W[:,0]
    flow["u"][:]   = W[:,1] / flow["rho"]
    flow["v"][:]   = W[:,2] / flow["rho"]
    flow["w"][:]   = W[:,3] / flow["rho"]
    flow["E"][:]   = W[:,4] / flow["rho"]

    if "cpg" in thermo:
        gamma = opts["gamma"]
        Pr = opts["Pr"]
        R = opts["R"]
        cp = gamma * R / (gamma - 1.)

        vel_mag_square = flow["u"]**2 + flow["v"]**2 + flow["w"]**2
        gamma_m1_rho = (gamma - 1.) * flow["rho"]

        flow["p"][:] = gamma_m1_rho * (flow["E"] - 0.5 * vel_mag_square)
        flow["E"][:] = flow["p"] / gamma_m1_rho + 0.5 * vel_mag_square
        flow["H"][:] = flow["E"] + flow["p"] / flow["rho"]
        flow["T"][:] = flow["p"] / (flow["rho"] * R)
        flow["c"][:] = sqrt(gamma * R * flow["T"])
        flow["mu"][:] = 1.45 * flow["T"]**1.5 * 1.e-6 / (flow["T"] + 110.)
        flow["cp"][:] = cp
        flow["k"][:] = cp * flow["mu"] / Pr
        flow["gamma"][:] = gamma
    else:
        print("error: unrecognized thermodynamics: " + thermo)
