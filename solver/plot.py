import matplotlib.pyplot as plt

def plot_residual_history(res_L2_norm_hist):
    for i in range(5):
        plt.semilogy(range(res_L2_norm_hist.shape[0]), res_L2_norm_hist[:,i] / res_L2_norm_hist[0,i])
    plt.legend(["continuity", "x-momentum", "y-momentum", "z-momentum", "energy"])
    plt.xlabel("Iterations")
    plt.ylabel("Normalized L2-norm of residuals")
    plt.show()
