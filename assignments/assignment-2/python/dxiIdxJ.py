import numpy as np

from psi import psi

def dxiIdxJ(elem, nodes_per_element, xi_1, xi_2, node_pos, elem_node):
    # 2-by-4 matrix containing global node coordinates of element
    elem_pos = np.zeros((2, nodes_per_element))
    # Matrix of grad(Psi_n) w.r.t. xi
    grad_mat = np.zeros((nodes_per_element, 2))
    for n in range(nodes_per_element):
        elem_pos[0, n] = node_pos[elem_node[elem, n], 0]
        elem_pos[1, n] = node_pos[elem_node[elem, n], 1]
        grad_mat[n, 0] = psi(n, 1, xi_1, xi_2)
        grad_mat[n, 1] = psi(n, 2, xi_1, xi_2)

    # Matrix of dxJ/dxiI
    jacobian_mat = np.matmul(elem_pos, grad_mat)
    jacobian = np.linalg.det(jacobian_mat)
    # Calculates dxi/dx which is the inverse of the jacobian matrix
    dxidx = np.linalg.inv(jacobian_mat)
    assert jacobian >= 0, "Jacobian must be positive"
    return jacobian, dxidx
