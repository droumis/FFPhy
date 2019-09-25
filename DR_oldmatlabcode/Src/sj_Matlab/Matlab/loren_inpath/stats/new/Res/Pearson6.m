% PEARSON6: Pearson Type VI curve
%
%     Usage: [xx,yy,mse] = pearson6(x,fx,mu,b1,b2,xmin,xmax)
%
%         x =       midpoints of histogram bars.
%         fx =      heights of histogram bars.
%         mu =      vector of moments.
%         b1,b2 =   beta values.
%         xmin =    lower extent of function on axis.
%         xmax =    upper extent of function on axis.
%         ----------------------------------------------------------
%         x,y =     coordinates of fitted function.
%         mse =     mean squared deviation between fx and predicted.
%

% RE Strauss, 3/4/01

function [xx,yy,mse] = pearson6(x,fx,mu,b1,b2,xmin,xmax)
  x = x - mu(1);

  r = 6*(b2-b1-1)/(6+3*b1-2*b2);
  rad = b1*(r+2)^2+16*(r+1);
  q1 = -((r-2)/2 - r*(r+2)/2*sqrt(b1/rad));
  q2 = (r-2)/2 + r*(r+2)/2*sqrt(b1/rad);
  a = abs(0.5*sqrt(mu(2))*sqrt(rad)) * sign(mu(3));
  A1 = (a*(q1-1))/((q1-1)-(q2-1));
  A2 = (a*(q2+1))/((q1-1)-(q2-1));

  y0 = 1;
  yy = real(y0 .* (1+x./A1).^(-q1) .* (1+x./A2).^q2);
  y0 = y0/sum(yy);
  yy = real(y0 .* (1+x./A1).^(-q1) .* (1+x./A2).^q2);
  mse = abs(mean((yy-fx).^2));

  xx = linspace(xmin-mu(1),xmax-mu(1))';
  yy = real(y0 .* (1+xx./A1).^(-q1) .* (1+xx./A2).^q2);

  xx = xx + mu(1);

  return;
