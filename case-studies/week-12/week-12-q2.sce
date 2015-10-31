clear;
stacksize('max');
//==============================================================================
// Creates sparse matrix using the compressed row storage format
//==============================================================================
function [val, col_ind, row_ptr] = MakeCSRMat(M)
  [num_rows, num_cols] = size(M);
  first_value = %T;
  for i = 1 : num_rows
    new_row = %T;
    for j = 1 : num_cols
      if first_value & M(i, j) ~= 0 then
        val = M(i, j);
        col_ind = j;
        row_ptr = 1;
        first_value = %F;
        new_row = %F;
      elseif M(i, j) ~= 0 then
        val($ + 1) = M(i, j);
        col_ind($ + 1) = j;
        if new_row then
          row_ptr($ + 1) = length(val);
          new_row = %F;
        end
      end
    end  // j
  end  // i
endfunction

//==============================================================================
// Performs matrix vector product of the form A * b = c, where A is a CSR matrix
//==============================================================================
function c = MatrixVectorProduct(A_val, A_col_ind, A_row_ptr, b)
  num_rows = length(A_row_ptr);
  num_vals = length(A_val);
  vec_size = length(b);
  c = zeros(num_rows, 1);
  row_ind = 0;
  for i = 1 : num_vals
    if find(i == A_row_ptr) then
      row_ind = row_ind + 1;
    end
    if A_col_ind(i) > vec_size then
      error('Matrix size too large');
    end
    c(row_ind) = c(row_ind) + A_val(i) * b(A_col_ind(i));
  end
endfunction

mat = [-1, 1, 0, 0;
       1, -2, 1, 0;
       0, 1, -2, 1;
       0, 0, 1, -1];
vec = [1; 2; 3; 4];

[csr_val, csr_col_ind, csr_row_ptr] = MakeCSRMat(mat);

c_vec = MatrixVectorProduct(csr_val, csr_col_ind, csr_row_ptr, vec);
