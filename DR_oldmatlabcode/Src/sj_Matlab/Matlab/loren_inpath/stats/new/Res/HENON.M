% Henon transform

function henon(x,y,nangles,npts)
  if (nargin < 3)
    nangles = 10;
  end;

  if (nargin < 4)
    npts = 1000;
  end;

  min_angle = 0;
  max_angle = pi/2;

  clf;
  hold on;

  for na = 1:nangles
    alpha = min_angle + (2*pi)*na/max_angle;
    pt = zeros(npts,2);
    pt(1,:) = [x y];

    for i = 2:npts
      xs = x;
      ys = x;
      d = ys - xs*xs;
       x = xs*cos(alpha) - d*sin(alpha);
       y = xs*sin(alpha) + d*cos(alpha);
       pt(i,:) = [x y];
    end;
min(pt)
max(pt)

    plot(pt(:,1),pt(:,2),'.');
  end;
%  axis([-1 1 -1 1]);
  hold off;

  return;
