clear;
stacksize('max');
//==============================================================================
// Creates a coordinate list sparse matrix. Appends non-zero values to the
// sparse matrix in the following format [row number, column number, value]
// Entries automatically sorted by row index and then column index.
// \param M dense matrix
// \return COO COO sparse matrix
//==============================================================================
function COO = MakeCOOMat(M)
  [num_rows, num_cols] = size(M);
  first_value = %T;
  for i = 1 : num_rows
    for j = 1 : num_cols
      if first_value & M(i, j) ~= 0 then
        COO(1, :) = [i, j, M(i, j)];
        first_value = %F;
      elseif M(i, j) ~= 0 then
        COO($ + 1, :) = [i, j, M(i, j)];
      end
    end  // j
  end  // i
endfunction

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
// Creates sparse matrix using the compressed row storage format. Inputs are
// symmetric matrices. Only stores upper half of matrix, ignores lower half.
//==============================================================================
function [val, col_ind, row_ptr] = MakeCSRMatSymmetric(M)
  [num_rows, num_cols] = size(M);
  first_value = %T;
  for i = 1 : num_rows
    new_row = %T;
    for j = i : num_cols
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
n = 4096;
path = pwd() + '\case-studies\week-12\';
A = read(path+'Amatrix', n, n);
B = read(path+'Bvector', n, 1);

mat = [-1, 1, 0, 0;
       1, -2, 1, 0;
       0, 1, -2, 1;
       0, 0, 1, -1];
vec = [2; 3; 4; 5];

coo_mat = MakeCOOMat(mat);
[csr_val, csr_col_ind, csr_row_ptr] = MakeCSRMat(mat);
[csrs_val, csrs_col_ind, csrs_row_ptr] = MakeCSRMatSymmetric(mat);

coo_A = MakeCOOMat(A);
[csr_val_A, csr_col_ind_A, csr_row_ptr_A] = MakeCSRMat(A);
[csrs_val_A, csrs_col_ind_A, csrs_row_ptr_A] = MakeCSRMatSymmetric(A);
