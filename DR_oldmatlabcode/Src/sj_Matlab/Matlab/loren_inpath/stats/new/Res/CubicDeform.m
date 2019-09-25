% CubicDeform: models a cubic-spline deformation in 1 dimension.
%
%     Usage: cubicdeform(basepts,targetpts,{ngrid})
%
%         basepts =  vector (length n) of landmark positions in 1D base form.
%         targetpts = matching vector of landmark positions in the target form.
%         ngrid =  optional number of grid points along line to be deformed [default = 20].
%

function cubicdeform(basepts,targetpts,ngrid)
  if (nargin < 3) ngrid = []; end;
  
  if (isempty(ngrid)) ngrid = 20; end;

  if (~isvector(basepts) | ~isvector(targetpts))
    error('  CubicDeform: base or target matrics must be vectors.');
  end;
  
  n = length(basepts);
  if (length(targetpts)~=n)
    error('  CubicDeform: base and target vectors must be of same length.');
  end;
  
  basemin = min(basepts) - 0.05*range(basepts);
  basemax = max(basepts) + 0.05*range(basepts);
  
  basegrid = linspace(basemin,basemax,ngrid);
  gridshift = zeros(1,ngrid);
  
  nshift = 100;
  base = linspace(basemin,basemax,nshift);
  shift = zeros(1,nshift);
  
  for i = 1:ngrid                           % Deform grid
    d = basegrid(i) - targetpts;
    gridshift(i) = sum(abs(d).^3);
  end;
  targetgrid = basegrid + gridshift;
  
  for i = 1:nshift                          % Deform "continuous" field
    d = base(i) - targetpts;
    shift = sum(abs(d).^3);
  end;
  target = base + shift;
  
  figure;
  plot(base,shift,'b',basepts,targetpts,'ro');
  putxlab('Position along form');
  putylab('Magnitude of shift');
  
  return;
  