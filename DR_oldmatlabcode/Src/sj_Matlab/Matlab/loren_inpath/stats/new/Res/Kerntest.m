% kerntest

clear;
close all;

load 'crds.dat';

if (1)
  n = 10;
  h = linspace(0.04,0.06,n);
  alpha = linspace(-1,-2,n);

  [H,A] = meshgrid(h,alpha);

  h = H(:);
  alpha = A(:);

  rl = zeros(size(h));
  for i = 1:length(h)
    rl(i) = kernlh([h(i),alpha(i)],crds);
  end;

  H = reshape(h,n,n);
  A = reshape(alpha,n,n);
  RL = reshape(rl,n,n);

  [c,ch] = contour(H,A,RL);
  clabel(c,ch);
  axis('square');
end;

if (0)
  n=20;
  x = linspace(0,1,n);
  y = linspace(0,1,n);

  [X,Y] = meshgrid(x,y);
  x = X(:);
  y = Y(:);

  h = mean(nndist(crds))
  fx = kernden([x y],crds,h);

  FX = reshape(fx,n,n);
  [c,ch] = contour(X,Y,FX);
  clabel(c,ch);
  axis('square');
end;


