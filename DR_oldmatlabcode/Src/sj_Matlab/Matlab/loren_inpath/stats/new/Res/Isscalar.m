% ISSCALAR: Returns 1 if the input matrix is a scalar, 0 otherwise.
%
%       b = isscalar(x)
%
%           x = input matrix.
%           -----------------
%           b = boolean flag.
%

% RE Strauss, 1/4/00

function b = isscalar(x)
  b = 1;

  if (isempty(x))
    b = 0;
  elseif (max(size(x)) > 1)
    b = 0;
  end;

  return;