% Fractal: fractal texture generator
%
%   Usage: fractal(target,fraction,iter)
%
%             target =   [n x 2] set of target points
%             fraction = fractional value between 0 and 1
%                           triangle:  1/2 or 2/3
%                           rectangle: 2/3
%                           pentagon:  2/3
%                           hexagon:   2/3
%             iter = number of iterations (plotted points)
%

function fractal(target,fraction,iter)
  [n,p] = size(target);

  plot(target(:,1),target(:,2),'ko');
  putbnd(target(:,1),target(:,2));
  hold on;

  p = mean(target);
  plot(p(1),p(2),'b.');

  for it = 1:iter
    r = ceil(rand*n);
    p = p + fraction*(target(r,:)-p);
    plot(p(1),p(2),'b.');
  end;

  return;

