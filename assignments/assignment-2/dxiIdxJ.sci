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

  // Create a 4-by-2 matrix containing global node coordinates of element
  elem_pos = zeros(nodesPerElement, 2);
  // Matrix of grad(dPsi_n / dxi_i)
  grad_mat = zeros(2, nodesPerElement);
  for n = 1 : nodesPerElement
    elem_pos(n, 1) = nodePos(elemNode(elem, n), 1);
    elem_pos(n, 2) = nodePos(elemNode(elem, n), 2);
    grad_mat(1, n) = psi(n, 1, xi1, xi2);
    grad_mat(2, n) = psi(n, 2, xi1, xi2);
  end  // n
  // matrix of dxJ/dxiI
  jacobian_mat = grad_mat * elem_pos;
  jacobian = det(jacobian_mat);
  dxidx = jacobian_mat \ grad_mat;
  // Error checking since jacobian must be positive
  if jacobian < 0 then
    error('Jacobian is negative.');
  end
endfunction
