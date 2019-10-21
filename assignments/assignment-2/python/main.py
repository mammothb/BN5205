import time

from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

from dxiIdxJ import dxiIdxJ 
from psi import psi

def main():
    start_time = time.time()
    # change resolution
    elems_per_side = 16

    side_length = 8
    nodes_per_side = elems_per_side + 1
    incr = side_length / elems_per_side
    nodes_per_element = 4
    num_nodes = nodes_per_side * nodes_per_side
    num_elems = elems_per_side * elems_per_side

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
    n = 0
    for j in range(nodes_per_side):
        for i in range(nodes_per_side):
            node_pos[n, 0] = i * incr
            node_pos[n, 1] = j * incr
            n += 1
    bc[:, 0] = 2

    ###################################################################
    # Generate elements
    ###################################################################
    e = 0
    for j in range(elems_per_side):
        for i in range(elems_per_side):
            elem_node[e, 0] = j * nodes_per_side + i
            elem_node[e, 1] = elem_node[e, 0] + 1
            elem_node[e, 2] = elem_node[e, 0] + nodes_per_side
            elem_node[e, 3] = elem_node[e, 0] + nodes_per_side + 1
            e += 1

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
    gauss_wei = np.array([0.25, 0.25, 0.25, 0.25])

    ###################################################################
    # Assemble global stiffness & mass matrices K & M
    ###################################################################
    K = np.zeros((num_nodes, num_nodes))
    M = np.zeros((num_nodes, num_nodes))
    f = np.zeros((num_nodes, 1))
    u = np.zeros((num_nodes, 1))

    for elem in range(num_elems):
        # initialize element stiffness and mass matrices
        EK = np.zeros((nodes_per_element, nodes_per_element))
        EM = np.zeros((nodes_per_element, nodes_per_element))
        ###############################################################
        # Create the element stiffness matrix
        ###############################################################
        for n in range(nodes_per_element):
            xi_1, xi_2 = gauss_pos[n]
            # Get jacobian and dxi/dx matrix
            J, dxidx = dxiIdxJ(elem, nodes_per_element, xi_1, xi_2, node_pos,
                               elem_node)
            # Matrix of grad(Psi_n) w.r.t. xi
            grad_psi = np.zeros((2, nodes_per_element))
            for nn in range(nodes_per_element):
                grad_psi[0, nn] = psi(nn, 1, xi_1, xi_2)
                grad_psi[1, nn] = psi(nn, 2, xi_1, xi_2)
            # Formula for B = dPsi/dxi * dxi/dx
            B = np.matmul(dxidx, grad_psi)
            # Setup Psi terms for Mass matrix
            psi_mat = np.array([psi(i, 0, xi_1, xi_2)
                                for i in range(nodes_per_element)])
            # Set up 4x4 matrix for Psi_n*Psi_m terms
            psi_big_mat = np.matmul(psi_mat.T, psi_mat)
            # Compute element stiffness matrix and element mass matric
            # with Gaussian quadrature
            EK += gauss_wei[n] * J * (D * np.matmul(B.T, B) - psi_big_mat)
            EM += gauss_wei[n] * J * psi_big_mat

        ###############################################################
        # Assemble EK and EM into global stiffness & mass matrices
        ###############################################################
        # Loop over local nodes of each element matrix
        for i in range(nodes_per_element):
            # Get global node number
            m = elem_node[elem, i]
            for j in range(nodes_per_element):
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
        np.arange(0, num_nodes, nodes_per_side),
        np.arange(nodes_per_side - 1, num_nodes, nodes_per_side),
        np.arange(1, nodes_per_side - 1),
        np.arange((nodes_per_side - 1) * nodes_per_side + 1, num_nodes - 1)
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

    middle = elems_per_side // 2 * nodes_per_side + elems_per_side // 2 + 1
    u[middle] = 1.0 * (8.0 / elems_per_side)**2
    print(time.time() - start_time)

    ###################################################################
    # Solve
    ###################################################################
    for __ in range(num_time_steps):
        u = np.linalg.solve(M + 0.5 * dt * K,
                            np.matmul(M - 0.5 * dt * K, u) +
                            dt * np.matmul(K, f))

    X, Y = np.meshgrid(node_pos[: nodes_per_side, 0],
                       node_pos[:: nodes_per_side, 1])
    soln = np.reshape(u, X.shape)
    fig = plt.figure()
    axes = fig.add_axes([0, 0, 1, 1], projection="3d")

    axes.plot_surface(X, Y, soln, cmap=cm.viridis)

    plt.show()

if __name__ == "__main__":
    main()
