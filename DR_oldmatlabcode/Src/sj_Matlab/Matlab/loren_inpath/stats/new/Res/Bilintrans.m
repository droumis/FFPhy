% BILINTRANS: Fits bilinear segmented regression line with a smooth logistic 
%             transition between line segments, using the 5-parameter model of 
%             Koops and Grossman (1993):
%
%               Y = a + b1*X - (b1-b2)*r*ln[1+exp((X-c)/r)]
%
%     Usage: [p,mse] = bilintrans(x,y,{doplot})
%
%         x =       vector of independent-variable values. 
%         y =       corresponding vector of dependent-variable values.
%         doplot =  optional boolean variable indicating that plot of data and 
%                     fitted line segments is to be produced [default = 0].
%         ---------------------------------------------------------------------
%         p =       parameters [a,b1,b2,c,r], where a is the 
%                     intercept, b_i are the slopes of the two linear segments, 
%                     c is the abscissa of the change-point, and r is the 
%                     transition-smoothness parameter (small r -- abrupt; 
%                     large r -- smooth).
%         mse =     mean residual sum-of-squares.
%

% Koops, W.J. & M. Grossman. 1993. Multiphasic allometry. Growth, Development & 
%   Aging 57:183-192.

% RE Strauss, 10/9/00

function [p,mse] = bilintrans(x,y,doplot)
  if (nargin < 3) doplot = []; end;

  if (isempty(doplot))
    doplot = 0;
  end;

  if (~isvector(x) | ~isvector(y))
    error(  'BILINTRANS: input must be vectors.');
  end;
  if (length(x) ~= length(y))
    error(  'BILINTRANS: input must be vectors of identical length.');
  end;

  x = x(:);                               % Convert to col vectors
  y = y(:);
  n = length(x);

  [x,i] = sort(x);                        % Sort data by abscissa values
  y = y(i);

  % Get initial parameter estimates from fitting a bilinear regression.

  [b,t] = bilinear(x,y);

  a0 = b(1);
  b10 = b(3);
  b20 = b(4);
  c0 = t;
  r0 = 0.1;

  p0 = [a0 b10 b20 c0 r0];                  % Initial parameter estimates
  p = fminsearch('bilintransf',p0,optimset('Display','off'),x,y);  % Optimization

  [mse,crds] = bilintransf(p,x,y);

  if (doplot)
    scatter(x,y);
    hold on;
    plot(crds(:,1),crds(:,2),'k');
    hold off;
  end;

  return;
