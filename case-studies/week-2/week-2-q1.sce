clear;
for i = 1:20
  if ~modulo(i, 15)
    printf("Fizz-Buzz\n");
  elseif ~modulo(i, 5)
    printf("Buzz\n");
  elseif ~modulo(i, 3)
    printf("Fizz\n");
  else
    printf("%d\n", i);
  end
end  // i
