% COMPLSET :  Given a list of variables, finds the observations having complete 
%             data (or vice versa).
%
%     Usage: obs = complset(X,indices,{getvars})
%
%         X =       [n x p] data matrix.
%         indices = vector of indices of variables.
%         getvars = optional boolean flag indicating, if true, that the list of 
%                     indices represents observations rather than variables 
%                     [default = 0].
%         ---------------------------------------------------------------------
%         obs =       column vector of corresponding observations (or variables, 
%                     if 'getvars' is true.
%

% RE Strauss, 12/18/00

function obs = complset(X,indices,getvars)
  if (nargin < 3) getvars = []; end;

  if (isempty(getvars))
    getvars = 0;
  end;

  indices = indices(:);

  if (getvars)
    X = X';
  end;
  [n,p] = size(X);

  if (max(indices) > p)
    error('  COMPLSET: maximum index greater than matrix size.');
  end;

  X = X(:,indices);                   % Subset matrix
  obs = find(isfinite(rowsum(X)));    % Complete subset of observatons

  return;