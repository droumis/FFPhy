% Gasket: randomly generate Sierpinski gasket within triangle
%   Kaplan & Glass, "Understanding Nonlinear Dynamics", p. 122.

% RE Strauss, 11/20/96
%   9/7/99 - changed plot colors for Matlab v5.

function gasket
  iter = 50000;

  x = zeros(iter,1);
  y = zeros(iter,1);

  Ax = 0; Ay = 0;
  Bx = 1; By = 0;
  Cx = cos(pi/3); Cy = sin(pi/3);
  plot([Ax Bx Cx Ax],[Ay By Cy Ay]);
  axis equal;
  axis off;
  hold on;

  disp('Digitize random point within triangle,');
  disp('  press <Enter>');
  [xx,yy] = ginput;

  x(1) = xx(1);
  y(1) = yy(1);

  r = ceil(3*rand(iter,1));       % Random integers [1,2,3]

  for i = 2:iter
    if (r(i)==1)
      x(i) = mean([x(i-1) Ax]);
      y(i) = mean([y(i-1) Ay]);
    elseif (r(i)==2)
      x(i) = mean([x(i-1) Bx]);
      y(i) = mean([y(i-1) By]);
    else
      x(i) = mean([x(i-1) Cx]);
      y(i) = mean([y(i-1) Cy]);
    end;
  end;

  plot(x,y,'k.');
  hold off;

  return;
