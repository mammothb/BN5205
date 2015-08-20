clear
function y = plotfun(x)
  y = sin(0.3 * x) + 0.2 * x
endfunction
x = [-10:0.1:36]
x2 = [-10:2:36]
plot(x, plotfun(x), 'r-')
plot(x2, plotfun(x2), 'k.')
