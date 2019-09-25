% LOGB: Logarithm of x to any base b.  
%
%     Usage: lx = logb(x,b)
%

% RE Strauss, 5/29/97

function lx = logb(x,b)
  lx = [];
  if (nargin < 2)
    error('LOGB: Invalid base');
  elseif (b<eps | abs(b-1)<eps)
    error('LOGB: Invalid base');
  else
    lx = log(x)/log(b);
  end;

  return;
