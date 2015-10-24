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
  for n = 1 : nodesPerElement
    elem_pos(n, 1) = nodePos(elemNode(elem, n), 1);
    elem_pos(n, 2) = nodePos(elemNode(elem, n), 2);
  end  // n
  // Matrix of grad(dPsi_n / dxi_i), -( ) are not expanded for clarity
  grad_mat = [-(1 - xi2), 1 - xi2,      -xi2, xi2;
              -(1 - xi1),    -xi1,   1 - xi1, xi1];
  // matrix of dxJ/dxiI
  jacobian_mat = grad_mat * elem_pos;

  jacobian = det(jacobian_mat);
  if jacobian < 0 then
    error('Jacobian is negative.');
  end
  dxidx = jacobian_mat \ grad_mat;
endfunction
