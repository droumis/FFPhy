% BILINEAR: Fits bilinear segmented regression line using reparameterized model 
%           of Chappell (1989).
%
%     Usage: [b,t,mse] = bilinear(x,y,{doplot})
%
%         x =       vector of independent-variable values. 
%         y =       corresponding vector of dependent-variable values.
%         doplot =  optional boolean variable indicating that plot of data and 
%                     fitted line segments is to be produced [default = 0].
%         ---------------------------------------------------------------------
%         b =       reparameterized parameters [a1,a2,b1,b2], where a_i are the 
%                     intercepts of the two linear segments and b_i are the 
%                     slopes.
%         t =       scalar indicating the abscissa of the change-point.
%         mse =     mean residual sum-of-squares.
%

% Chappell, R. 1989. Fitting bent lines to data, with applications to 
%   allometry.  J. Theor. Biol. 138:235-256.
%   [NB error on p 239: min(0,x-t) should be max(0,x-t)].

% RE Strauss, 10/8/00

function [b,t,mse] = bilinear(x,y,doplot)
  if (nargin < 3) doplot = []; end;

  if (isempty(doplot))
    doplot = 0;
  end;

  if (~isvector(x) | ~isvector(y))
    error('  BILINEAR: input must be vectors.');
  end;
  if (length(x) ~= length(y))
    error('  BILINEAR: input must be vectors of identical length.');
  end;

  x = x(:);                               % Convert to col vectors
  y = y(:);
  n = length(x);

  sx = sort(x);

  t = fminbnd('bilinearf',[sx(2)+eps],[sx(n-1)-eps],optimset('Display','off'),x,y);
  [mse,b,crds] = bilinearf(t,x,y);

  if (doplot)
    scatter(x,y);
    hold on;
    plot(crds(:,1),crds(:,2),'k');
    hold off;
  end;

  return;
