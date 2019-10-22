import numpy as np

from psi import psi

def dxiIdxJ(elem, node_per_elem, xi_1, xi_2, node_pos, elem_node):
    # 2-by-4 matrix containing global node coordinates of element
    nodes = [elem_node[elem, num] for num in range(node_per_elem)]
    elem_pos = np.array([[node_pos[node, x] for node in nodes]
                         for x in range(2)])
    # Matrix of grad(Psi_n) w.r.t. xi
    grad_mat = np.array([[psi(num, der, xi_1, xi_2) for der in [1, 2]]
                         for num in range(node_per_elem)])
    # Matrix of dxJ/dxiI
    jacobian_mat = elem_pos @ grad_mat
    jacobian = np.linalg.det(jacobian_mat)
    # Calculates dxi/dx which is the inverse of the jacobian matrix
    dxidx = np.linalg.inv(jacobian_mat)
    assert jacobian >= 0, "Jacobian must be positive"
    return jacobian, dxidx
