import itertools
import time

from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

from dxiIdxJ import dxiIdxJ 
from psi import psi

def main():
    # start_time = time.time()
    # change resolution
    elem_per_side = 16

    side_length = 8
    node_per_side = elem_per_side + 1
    incr = side_length / elem_per_side
    node_per_elem = 4
    num_nodes = node_per_side * node_per_side
    num_elems = elem_per_side * elem_per_side

    node_pos = np.zeros((num_nodes, 2))
    elem_node = np.zeros((num_elems, 4), dtype=np.int32)
    bc = np.zeros((num_nodes, 2))

    ###################################################################
    # Set up model / simulation parameters
    ###################################################################
    D = 0.1
    sim_time = 10.0  # [s]
    dt = 0.05  # [s]
    num_time_steps = round(sim_time / dt)

    ###################################################################
    # Generate nodes & boundary conditions
    ###################################################################
    x_pos = np.arange(node_per_side) * incr
    for j in range(node_per_side):
        n_start = j * node_per_side
        n_stop = (j + 1) * node_per_side
        node_pos[n_start : n_stop, 0] = x_pos
        node_pos[n_start : n_stop, 1] = j * incr
    bc[:, 0] = 2

    ###################################################################
    # Generate elements
    ###################################################################
    row = np.arange(elem_per_side)
    for j in range(elem_per_side):
        e_start = j * elem_per_side
        e_stop = (j + 1) * elem_per_side
        elem_node[e_start : e_stop, 0] = j * node_per_side + row
    elem_node[:, 1] = elem_node[:, 0] + 1
    elem_node[:, 2] = elem_node[:, 0] + node_per_side
    elem_node[:, 3] = elem_node[:, 0] + node_per_side + 1

    ###################################################################
    # Set up numerical integration points and weights
    ###################################################################
    num_gauss_points = 4
    gp_1 = 0.5 - 1 / (2 * 3**0.5)
    gp_2 = 0.5 + 1 / (2 * 3**0.5)

    gauss_pos = np.array([[gp_1, gp_1],
                          [gp_2, gp_2],
                          [gp_1, gp_2],
                          [gp_2, gp_2]])
    gauss_wei = np.array([0.25] * num_gauss_points)

    ###################################################################
    # Assemble global stiffness & mass matrices K & M
    ###################################################################
    K = np.zeros((num_nodes, num_nodes))
    M = np.zeros((num_nodes, num_nodes))
    f = np.zeros((num_nodes, 1))
    u = np.zeros((num_nodes, 1))

    for elem in range(num_elems):
        # initialize element stiffness and mass matrices
        EK = np.zeros((node_per_elem, node_per_elem))
        EM = np.zeros((node_per_elem, node_per_elem))
        ###############################################################
        # Create the element stiffness matrix
        ###############################################################
        for n in range(node_per_elem):
            xi_1, xi_2 = gauss_pos[n]
            # Get jacobian and dxi/dx matrix
            J, dxidx = dxiIdxJ(elem, node_per_elem, xi_1, xi_2, node_pos,
                               elem_node)
            # Matrix of grad(Psi_n) w.r.t. xi
            grad_psi = np.array([[psi(num, der, xi_1, xi_2)
                                  for num in range(node_per_elem)]
                                 for der in [1, 2]])
            # Formula for B = dPsi/dxi * dxi/dx
            B = dxidx @ grad_psi
            # Setup Psi terms for Mass matrix
            psi_mat = np.array([psi(num, 0, xi_1, xi_2)
                                for num in range(node_per_elem)])
            # Set up 4x4 matrix for Psi_n*Psi_m terms
            psi_big_mat = psi_mat.T @ psi_mat
            # Compute element stiffness matrix and element mass matric
            # with Gaussian quadrature
            EK += gauss_wei[n] * J * (D * B.T @ B - psi_big_mat)
            EM += gauss_wei[n] * J * psi_big_mat

        ###############################################################
        # Assemble EK and EM into global stiffness & mass matrices
        ###############################################################
        # Loop over local nodes of each element matrix
        for i in range(node_per_elem):
            # Get global node number
            m = elem_node[elem, i]
            for j in range(node_per_elem):
                # Get the other global node number
                n = elem_node[elem, j]
                # Assembles element matrix into global matrix, sum up
                # values when global matrix nodes coincide
                K[m, n] += EK[i, j]
                M[m, n] += EM[i, j]

    ###################################################################
    # Apply boundary & initial conditions
    ###################################################################
    ext_nodes = np.concatenate((
        np.arange(0, num_nodes, node_per_side),
        np.arange(node_per_side - 1, num_nodes, node_per_side),
        np.arange(1, node_per_side - 1),
        np.arange((node_per_side - 1) * node_per_side + 1, num_nodes - 1)
        ))
    for n in range(num_nodes):
        # Essential boundary
        if bc[n, 0] == 1:
            K[n, :] = 0
            K[n, n] = 1
            f[n] = bc[n, 1]
        # Natural boundary (integrated fluxes only)
        # only apply BC to RHS if node is at the edge / an external
        # node since flux for internal nodes is zero
        elif n in ext_nodes and bc[n, 0] == 2:
            f[n] = bc[n, 1]

    middle = elem_per_side // 2 * node_per_side + elem_per_side // 2 + 1
    u[middle] = 1.0 * (8.0 / elem_per_side)**2

    ###################################################################
    # Solve
    ###################################################################
    for __ in range(num_time_steps):
        u = np.linalg.solve(M + dt / 2 * K, (M - dt / 2 * K) @ u + dt * K @ f)

    X, Y = np.meshgrid(node_pos[: node_per_side, 0],
                       node_pos[:: node_per_side, 1])
    soln = np.reshape(u, X.shape)
    fig = plt.figure()
    axes = fig.add_axes([0, 0, 1, 1], projection="3d")

    axes.plot_surface(X, Y, soln, cmap=cm.viridis)

    plt.show()

if __name__ == "__main__":
    main()
