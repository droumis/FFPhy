% SETRANGE: Alters the ranges of a matrix, by column, given min and max values.
%
%     Usage: Y = setrange(X,{minvals},{maxvals})
%
%           X =       [r x c] matrix.
%           minvals = optional vector (length c) of new minimum values for each 
%                       column [default = zeros].
%           maxvals = optional vector (length c) of new maximum values for each 
%                       column [default = ones].
%           -------------------------------------------------------------------
%           Y = transformed matrix.
%

% RE Strauss, 6/7/00

function Y = setrange(X,minvals,maxvals)
  if (nargin < 2) minvals = []; end;
  if (nargin < 3) maxvals = []; end;

  [r,c] = size(X);
  if (isempty(minvals))
    minvals = zeros(1,c);
  end;
  if (isempty(maxvals))
    maxvals = ones(1,c);
  end;

  if (length(minvals)~=c | length(maxvals)~=c)
    error('  SETRANGE: lengths of min & max vectors not compatible with X.');
  end;

  minvals = ones(r,1)*(minvals(:)');
  maxvals = ones(r,1)*(maxvals(:)');
  rngvals = maxvals - minvals;

  curmin = ones(r,1)*min(X);
  curmax = ones(r,1)*max(X);
  currng = curmax - curmin;

  Y = (X-curmin)./currng;
  Y = Y.*rngvals + minvals;

  return;
