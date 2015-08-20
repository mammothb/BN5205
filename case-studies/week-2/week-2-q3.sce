clear;
LHS = [3 5 1 0;
       1 0 1 0;
       1 3 0 1;
       -4 0 6 1];
RHS = [6; 63; 20; 3];
result = LHS \ RHS;
for i = result
  printf("%d ", i);
end  // i
