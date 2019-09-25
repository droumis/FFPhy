% QUADRATIC:  Fits quadratic regression line.
%
%     Usage: [b,mse] = quadratic(x,y,{doplot})
%
%         x =       vector of independent-variable values. 
%         y =       corresponding vector of dependent-variable values.
%         doplot =  optional boolean variable indicating that plot of data and 
%                     fitted line segments is to be produced [default = 0].
%         ---------------------------------------------------------------------
%         b =       parameters [b0,b1,b2].
%         mse =     mean residual sum-of-squares.
%

% RE Strauss, 10/9/00

function [b,mse] = quadratic(x,y,doplot)
  if (nargin < 3) doplot = []; end;

  if (isempty(doplot))
    doplot = 0;
  end;

  if (~isvector(x) | ~isvector(y))
    error(  'QUADRATIC: input must be vectors.');
  end;
  if (length(x) ~= length(y))
    error(  'QUADRATIC: input must be vectors of identical length.');
  end;

  x = x(:);                               % Convert to col vectors
  y = y(:);
  n = length(x);

  X = [ones(n,1) x x.*x];
  b = inv(X'*X) * X'*y;                   % Multiple regression

  e = y - (b(1) + b(2)*x + b(3)*x.*x);
  mse = (e'*e)/(n-3);

  if (doplot)
    xx = linspace(min(x),max(x));
    yy = b(1) + b(2)*xx + b(3)*xx.*xx;

    scatter(x,y);
    hold on;
    plot(xx,yy,'k');
    hold off;
  end;

  return;


