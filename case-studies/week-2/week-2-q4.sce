clear
nums = [1, 1];
ans = 0;
while sum(nums) < 4e6
  if ~modulo(sum(nums), 2)
    ans = ans + sum(nums);
  end
  nums = [nums(2), sum(nums)];
end
printf("%d", ans)
