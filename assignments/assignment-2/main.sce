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
clf();

// Path to your sci files - adjust as needed
// path = 'C:\Users\biebml\Documents\';
path = pwd();

// Load necessary functions into scilab
exec(path+'\assignments\assignment-2\psi.sci');
exec(path+'\assignments\assignment-2\dxiIdxJ.sci');

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
elemsPerSide = 64;

sideLength = 8;
nodesPerSide = elemsPerSide+1;
incr = sideLength/elemsPerSide;
nodesPerElement = 4;
numNodes = nodesPerSide * nodesPerSide;
numElems = elemsPerSide * elemsPerSide;
nodePos(1:numNodes,2) = 0;
elemNode(1:numElems,4) = 0;
BCs(1:numNodes,2) = 0;

//==============================================================================
// Set up model / simulation parameters
//==============================================================================
D = 0.1;  // Diffusivity of domain
sim_time = 10.0;  // [s], total simulation time
dt = 0.05;  // [s], time step
num_time_steps = round(sim_time / dt);  // total number of time steps
// Warns if sim_time is not completely divisible by dt, either due to incorrect
// choice of dt or due to precision error. Example when sim_time = 0.3 and
// dt = 0.1
if modulo(sim_time, dt) > 1e-20 then
  warning('Inappropriate time step size');
end

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
  for n = 1 : nodesPerElement
    // Get values for Gauss points
    xi1 = gaussPos(n, 1);
    xi2 = gaussPos(n, 2);
    // Calculate jacobian and dxi/dx matrix
    [J, dxidx] = dxiIdxJ(elem, nodesPerElement, xi1, xi2, nodePos, elemNode);
    // Set up Psi terms for Mass matrix
    psi_mat = [psi(1, 0, xi1, xi2), psi(2, 0, xi1, xi2), psi(3, 0, xi1, xi2),...
        psi(4, 0, xi1, xi2)];
    // Set up 4x4 matrix for Psi_n*Psi_m terms
    psi_big_mat = psi_mat' * psi_mat;
    // Compute element stiffness matrix and element mass matrix with Gaussian
    // quadrature
    EK = EK + gaussWei(n) * J .* (D .* dxidx' * dxidx - psi_big_mat);
    EM = EM + gaussWei(n) * J .* psi_big_mat;
  end  // nn

  //===========================================================================
  // Assemble EK & EM into the global stiffness & mass matrices
  //===========================================================================

  //*** Q. Write code to calculate assemble the element stiffness
  //***    matrix you have just created into the global stiffness matrix here.

  // Loop over local nodes of each element matrix
  for i = 1 : nodesPerElement
    m = elemNode(elem, i);  // Get the global node number
    for j = 1 : nodesPerElement
      n = elemNode(elem, j);  // Get the other global node number
      // printf('%d, %d\n', m, n);\
      // Assemles element matrix into global matrix, sum up values when global
      // matrix nodes coincide.
      K(m, n) = K(m, n) + EK(i, j);
      M(m, n) = M(m, n) + EM(i, j);
    end  // j
  end  // i
end //elem

//=============================================================================
// Apply boundary & initial conditions
//=============================================================================

//*** Q. Write code to apply the boundary and initial conditions to the problem
//***    here.
//***    NOTE: the overwriting method should be used for essential BCs.

// Vector containing global node number of external/edge nodes in the following
// format: ext_nodes = [left edge, right edge, bottom edge, top edge]
ext_nodes = [1 : nodesPerSide : numNodes,...
             nodesPerSide : nodesPerSide : numNodes,...
             2 : nodesPerSide - 1,...
             (nodesPerSide - 1) * nodesPerSide + 2 : numNodes - 1];
// Loop over all nodes
for n = 1 : numNodes
  // Essential boundary
  if BCs(n, 1) == 1 then
    // Overwrite row with zeros and then put 1 at the position with essential BC
    K(n, :) = zeros(1, numNodes);
    K(n, n) = 1;
    // Set BC value for RHS vector
    f(n) = BCs(n, 2);
  // Natural boundary (integrated fluxes only)
  // Only apply BC to RHS if node is at the edge / an external node since flux
  // for internal nodes is zero
  elseif find(ext_nodes == n) & BCs(n, 1) == 2 then
    f(n) = BCs(n, 2);
  end
end  // n
// Initial condition of u(4, 4, t = 0) = 1. Only works for odd nodesPerSide
// there will not be a single center node when nodesPerSide is even.
u(elemsPerSide / 2 * nodesPerSide + elemsPerSide / 2 + 1) = 1.0;

//=============================================================================
// Solve
//=============================================================================

//*** Q. Write code to solve the resulting matrix system over time here.

// Solves using Crank-Nicolson-Galerkin method. Equation has the following form:
// [M + \theta * \Delta t * K] * u^{n + 1} =
//     [M - (1 - \theta) * \Delta t * K] * u^n + \Delta t * K * f
// where \theta = 0.5
// 2nd order accurate in time domain and unconditionably stable since it's
// implicit in time.
// Able to obtain stable solution when dt = 0.5 whereas Forward-Euler shows
// instability
for t = 1 : num_time_steps
  // Only one solution vector since we are only interested in the end result.
  u = (M + 0.5 * dt .* K) \ ((M - 0.5 * dt .* K) * u + dt .* K * f);
end  // t

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
surf(soln);
