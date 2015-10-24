//======================================================
// Function to calculate basis functions and derivatives
//
// num = local node number
// der = derivative: 0 = function; 1 = d/dxi_1; 2 = d/dxi_2
// xi1 = xi1 coordinate (range 0-1)
// xi2 = xi2 coordinate (range 0-1)
//======================================================

function result = psi (num, der, xi1, xi2 )
    
  if( der == 0 )
    select num 
      case 1 then
        result = (1-xi1) * (1-xi2);
      case 2 then
        result = xi1 * (1-xi2);
      case 3 then
        result = (1-xi1) * xi2;
      case 4 then
        result = xi1 * xi2;
      else
        result = 0;
    end
  elseif( der == 1 )
    select num 
      case 1 then
        result = -1 * (1-xi2);
      case 2 then
        result = (1-xi2);
      case 3 then
        result = -xi2;
      case 4 then
        result = xi2;
      else
        result = 0;
    end
  elseif( der == 2 )
    select num 
      case 1 then
        result = -1 * (1-xi1);
      case 2 then
        result = -xi1;
      case 3 then
        result = (1-xi1);
      case 4 then
        result = xi1;
      else
        result = 0;
    end
  else
    result = 0;
  end  // der

endfunction
