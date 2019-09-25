% CAUCHYRND: Returns random samples from the 2-parameter Cauchy distribution 
%            C(a,b).
%
%     Usage: rnd = Cauchyrnd(a,b,{[r,c]})
%
%         a =     location parameters (matrix or scalar), range -Inf to +Inf.
%         b =     scale parameters (matrix or scalar), range > 0.
%         [r,c] = optional vector of rows and columns for resulting matrix.  If 
%                   not supplied, then the size of the resulting matrix is the 
%                   common size of the parameter matrices which, if not both 
%                   scalars, must be of identical size.
%         ---------------------------------------------------------------------
%         rnd =   matrix of random values.
%

% RE Strauss, 6/7/00

function rnd = Cauchyrnd(location,scale,matsize)
  if (nargin < 3) matsize = []; end;

  if (~isempty(matsize))
    [isvect,len] = isvector(matsize);
    if (~isvect | len~=2)
      error('  CAUCHYRND: invalid matrix-size vector.');
    end;
    r = matsize(1);
    c = matsize(2);
  end;

  [samesize,location,scale] = commonsize(location,scale,zeros(r,c));
  if (~samesize)
    error('  CAUCHYRND: incompatible parameter matrices.');
  end;
  
  location = location(:)';
  scale = scale(:)';
  
  rnd = NaN*ones(size(location));
  k = find(isfinite(location) & isfinite(scale) & (scale>0));
  if (~isempty(k))
    rnd(k) = location(k) - cot(pi*rand(1,length(k))).*scale(k);
  end;

  rnd = reshape(rnd,r,c);
  
  return;

