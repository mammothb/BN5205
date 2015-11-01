//=============================================================================
// Function to calculate dxi_i/dx_j and the element Jacobian J
//
// elem = global element number
// nodesPerElement = number of nodes per element
// xi1 = xi1 coordinate (range 0-1)
// xi2 = xi2 coordinate (range 0-1)
// nodePos(numNodes, dim) = matrix of node coordinates
// elemNode = list of nodes for each element
//=============================================================================

function [jacobian, dxidx] = dxiIdxJ ( elem, nodesPerElement, xi1, xi2, nodePos, elemNode )

  //*** Q. Write code to calculate the term: dxi_i/dx_j
  //***    and the transformation Jacobian for an
  //***    arbitrary element in 2D space here.

  // Create a 2-by-4 matrix containing global node coordinates of element
  elem_pos = zeros(2, nodesPerElement);
  // Matrix of grad(Psi_n) w.r.t. xi
  grad_mat = zeros(nodesPerElement, 2);
  for n = 1 : nodesPerElement
    elem_pos(1, n) = nodePos(elemNode(elem, n), 1);
    elem_pos(2, n) = nodePos(elemNode(elem, n), 2);
    grad_mat(n, 1) = psi(n, 1, xi1, xi2);
    grad_mat(n, 2) = psi(n, 2, xi1, xi2);
  end  // n
  // Matrix of dxJ/dxiI
  jacobian_mat = elem_pos * grad_mat;
  jacobian = det(jacobian_mat);
  // Calculates dPsi_n/dxi_i * dxi_i/dx_j = inv(jacobian_mat) * grad(Psi_n) and
  // output that as dxi_i/dxi_j so we don't have to loop through element local
  // nodes to construct grad(Psi_n) again
  dxidx = jacobian_mat \ grad_mat';
  // Error checking since jacobian must be positive
  if jacobian < 0 then, error('Jacobian is negative.'); end
endfunction
