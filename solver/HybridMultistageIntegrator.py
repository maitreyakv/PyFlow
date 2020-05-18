from numpy import sqrt, abs

from solver.Integrator import Integrator

# TODO: Maybe make these options instead of hardcoded
# Five stage coefficients for Central Scheme
alpha_1 = 0.25
alpha_2 = 0.1667
alpha_3 = 0.375
alpha_4 = 0.5
alpha_5 = 1.
beta_3 = 0.56
beta_5 = 0.44

# CFL number for central scheme
CFL = 3.6

# Spectral radii constant for central scheme
C = 4.

# TODO: Add doc for class
class HybridMultistageIntegrator(Integrator):

    # Constructor for HybridMultistageIntegrator
    def __init__(self, interior_cells):
        # Superclass Constructor
        super().__init__()

        # Allocate maps for stage variables and residuals
        self.W_0 = {cell: None for cell in interior_cells}
        self.Rc_0 = {cell: None for cell in interior_cells}
        self.Rd_0 = {cell: None for cell in interior_cells}
        self.W_1 = {cell: None for cell in interior_cells}
        self.Rc_1 = {cell: None for cell in interior_cells}
        self.W_2 = {cell: None for cell in interior_cells}
        self.Rc_2 = {cell: None for cell in interior_cells}
        self.Rd_2_0 = {cell: None for cell in interior_cells}
        self.W_3 = {cell: None for cell in interior_cells}
        self.Rc_3 = {cell: None for cell in interior_cells}
        self.W_4 = {cell: None for cell in interior_cells}
        self.Rc_4 = {cell: None for cell in interior_cells}
        self.Rd_4_2 = {cell: None for cell in interior_cells}

        # Compute projections of control volumes onto principle planes
        self.area_proj = {cell: 0.5 * sum([abs(face.area * face.n_) for face in cell.faces]) for cell in interior_cells}

    # Implements integration using hybrid multi-stage explicit (Runge-Kutta)
    def integrate(self, sim, t, dt=None):
        # Compute initial residuals
        sim.compute_residuals(t)

        # Compute time step
        dt = 1e10
        for cell in sim.interior_cells:
            spec_rad_c = (abs(cell.flow.v_) + cell.flow.c) * self.area_proj[cell]
            spec_rad_v = max(4. / (3. * cell.flow.rho), cell.flow.gamma / cell.flow.rho) * (cell.flow.k / cell.flow.cp) * self.area_proj[cell]**2 / cell.vol
            dt_local = CFL * cell.vol / (sum(spec_rad_c) + C * sum(spec_rad_v))
            dt = min(dt, dt_local)

        # Perform first stage
        for cell in sim.interior_cells:
            self.Rc_0[cell] = cell.Rc
            self.Rd_0[cell] = cell.Rd
            self.W_0[cell] = cell.flow.W_
            self.W_1[cell] = self.W_0[cell] - alpha_1 * dt * (self.Rc_0[cell] - self.Rd_0[cell]) / cell.vol

        # Perform second stage
        for cell in sim.interior_cells:
            cell.flow.W_ = self.W_1[cell]
        sim.compute_residuals(t)
        for cell in sim.interior_cells:
            self.Rc_1[cell] = cell.Rc
            self.W_2[cell] = self.W_0[cell] - alpha_2 * dt * (self.Rc_1[cell] - self.Rd_0[cell]) / cell.vol

        # Perform third stage
        for cell in sim.interior_cells:
            cell.flow.W_ = self.W_2[cell]
        sim.compute_residuals(t)
        for cell in sim.interior_cells:
            self.Rc_2[cell] = cell.Rc
            self.Rd_2_0[cell] = beta_3 * cell.Rd + (1. - beta_3) * self.Rd_0[cell]
            self.W_3[cell] = self.W_0[cell] - alpha_3 * dt * (self.Rc_2[cell] - self.Rd_2_0[cell]) / cell.vol

        # Perform fourth stage
        for cell in sim.interior_cells:
            cell.flow.W_ = self.W_3[cell]
        sim.compute_residuals(t)
        for cell in sim.interior_cells:
            self.Rc_3[cell] = cell.Rc
            self.W_4[cell] = self.W_0[cell] - alpha_4 * dt * (self.Rc_3[cell] - self.Rd_2_0[cell]) / cell.vol

        # Perform fifth stage
        for cell in sim.interior_cells:
            cell.flow.W_ = self.W_4[cell]
        sim.compute_residuals(t)
        for cell in sim.interior_cells:
            self.Rc_4[cell] = cell.Rc
            self.Rd_4_2[cell] = beta_5 * cell.Rd + (1. - beta_5) * self.Rd_2_0[cell]
            cell.flow.W_ = self.W_0[cell] - alpha_5 * dt * (self.Rc_4[cell] - self.Rd_4_2[cell]) / cell.vol

        # Compute final residual L2 norm
        res_L2_norm = sqrt(sum([(self.Rc_4[cell] - self.Rd_4_2[cell])**2 for cell in sim.interior_cells]))

        return dt, res_L2_norm
