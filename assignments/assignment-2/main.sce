//=============================================================================
// Program to solve 2D Reaction Diffusion equation using Galerkin FE
//
// by   :
// date :
//=============================================================================

// Clear all existing variables
clear;

// Do nothing if a function is redefined
funcprot(0);

// This is needed to solve higher resolution meshes
stacksize('max');

// Clear previous figures, if any
// clf();

// Path to your sci files - adjust as needed
// path = 'C:\Users\biebml\Documents\';
path = pwd();

// Load necessary functions into scilab
exec(path+'\psi.sci');
exec(path+'\dxiIdxJ.sci');

//=============================================================================
// Set up the mesh
//
// numNodes = total number of global nodes
// numElems = total number of global elements
//
// nodePos(n,1) = x position of node n
// nodePos(n,2) = y position of node n
//
// elemNode(e,1) = first global node in element e (local node 1)
// elemNode(e,2) = second global node in element e (local node 2)
// elemNode(e,3) = third global node in element e (local node 3)
// elemNode(e,4) = fourth global node in element e (local node 4)
//
// BC(n,1) = type of boundary condition at node n
//           BC(n,1) = 1 => essential (potential)
//           BC(n,1) = 2 => natural (flux)
// BC(n,2) = value of the boundary condition at node n
//=============================================================================

// Change resolution using this variable
elemsPerSide = 8;

sideLength = 8;
nodesPerSide = elemsPerSide+1;
incr = sideLength/elemsPerSide;
nodesPerElement = 4;
numNodes = nodesPerSide * nodesPerSide;
numElems = elemsPerSide * elemsPerSide;
nodePos(1:numNodes,2) = 0;
elemNode(1:numElems,4) = 0;
BCs(1:numNodes,2) = 0;

//=============================================================================
// Generate nodes & boundary conditions
//=============================================================================
n = 1;
for j = 1:nodesPerSide
  for i = 1:nodesPerSide
    nodePos(n,1) = (i-1)*incr;
    nodePos(n,2) = (j-1)*incr;
    BCs(n,1) = 2;
    BCs(n,2) = 0;
    n = n + 1;
  end //i
end //j

//=============================================================================
// Generate elements
//=============================================================================
e = 1;
for j = 1:elemsPerSide
  for i = 1:elemsPerSide
    elemNode(e,1) = (j-1)*nodesPerSide+i;
    elemNode(e,2) = elemNode(e,1)+1;
    elemNode(e,3) = elemNode(e,1)+nodesPerSide;
    elemNode(e,4) = elemNode(e,1)+nodesPerSide+1;
    e = e + 1;
  end //i
end //j

//=============================================================================
// Set up the numerical integration points 'gaussPos' and weights 'gaussWei'
//=============================================================================
numGaussPoints = 4; // 2x2

gp1 = 0.5 - ( 1 / ( 2 * sqrt( 3 )));
gp2 = 0.5 + ( 1 / ( 2 * sqrt( 3 )));

gaussPos = [ gp1 gp1;
             gp2 gp1;
             gp1 gp2;
             gp2 gp2 ];

gaussWei = [ 0.25;
             0.25;
             0.25;
             0.25 ];
//=============================================================================
// Assemble the global stiffness & mass matrices K & M
//=============================================================================
K = zeros(numNodes,numNodes);
M = zeros(numNodes,numNodes);
f = zeros(numNodes,1);
u = zeros(numNodes,1);
D = 0.1;  // Diffusivity of domain

// Loop over elements
for elem = 1:numElems

  // Initialise element stiffness (EK) & mass (EM)
  EK = zeros(nodesPerElement,nodesPerElement);
  EM = zeros(nodesPerElement,nodesPerElement);

  //===========================================================================
  // Create the element stiffness matrix
  //===========================================================================

  //*** Q. Write code to calculate the element stiffness & mass
  //***    matrices for an arbitrary element in 2D space here.
  // Loop over local nodes of each element
  for nn = 1 : nodesPerElement
    // Get values for Gauss points
    xi1 = gaussPos(nn, 1);
    xi2 = gaussPos(nn, 2);
    // Calculate jacobian and dxi/dx matrix
    [J, dxidx] = dxiIdxJ(elem, nodesPerElement, xi1, xi2, nodePos, elemNode);
    // Compute element stiffness matrix with Gaussian quadrature
    EK = EK + gaussWei(nn) * D * J .* dxidx' * dxidx;
  end
  disp(EK);

  //===========================================================================
  // Assemble EK & EM into the global stiffness & mass matrices
  //===========================================================================

  //*** Q. Write code to calculate assemble the element stiffness
  //***    matrix you have just created into the global stiffness matrix here.

end //elem

//=============================================================================
// Apply boundary & initial conditions
//=============================================================================

//*** Q. Write code to apply the boundary and initial conditions to the problem
//***    here.
//***    NOTE: the overwriting method should be used for essential BCs.

//=============================================================================
// Solve
//=============================================================================

//*** Q. Write code to solve the resulting matrix system over time here.

//=============================================================================
// Surface plot of the solution
//=============================================================================
n = 1;
for j = 1:nodesPerSide
  for i = 1:nodesPerSide
    soln(i,j) = u(n);
    n = n + 1;
  end //i
end //j
// surf(soln);
