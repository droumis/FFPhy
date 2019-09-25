% PEARSON1: Pearson Type I curve
%
%     Usage: [xx,yy,mse] = pearson1(x,fx,mu,b1,b2)
%
%         x =       midpoints of histogram bars.
%         fx =      heights of histogram bars.
%         mu =      vector of moments.
%         b1,b2 =   beta values.
%         ----------------------------------------------------------
%         x,y =     coordinates of fitted function.
%         mse =     mean squared deviation between fx and predicted.
%

% RE Strauss, 3/4/01

function [xx,yy,mse] = pearson1(x,fx,mu,b1,b2)

  r = 6*(b2-b1-1)/(6+3*b1-2*b2);
  mode = mu(1) - 0.5*(mu(3)/mu(2))*((r+2)/(r-2));
  x = x - mode;

  rad = sqrt(b1/(b1*(r+2)^2+16*(r+1)));
  m1 = 0.5*((r-2) - r*(r+2)*rad);
  m2 = 0.5*((r-2) + r*(r+2)*rad);
  asum = 0.5*sqrt(mu(2))*sqrt(b1*(r+2)^2+16*(r+1));
  a1 = asum/(1+m2/m1);
  a2 = asum-a1;

  y0 = 1;
  yy = y0*((1+x/a1).^m1).*((1-x/a2).^m2);
  y0 = y0/sum(yy);
  yy = y0*((1+x/a1).^m1).*((1-x/a2).^m2);
  mse = abs(mean((yy-fx).^2));

  xx = linspace(-a1,a2)';
  yy = y0*((1+xx/a1).^m1).*((1-xx/a2).^m2);

  xx = xx+mode;

  return;
