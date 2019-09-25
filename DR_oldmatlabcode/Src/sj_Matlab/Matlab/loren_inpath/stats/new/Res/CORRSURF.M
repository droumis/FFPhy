% CORRSURF: Find the fitness surface of a correlation matrix as defined by 
%           two particular orthogonal components of the full permutation
%           space.
%
%     Usage: corrsurf(r)
%
%           r = square symmetric matrix
%

% RE Strauss, 4/3/95

function corrsurf(r)
  maxiter = 10;
  N = size(r,1);
  w = corradj(N);                     % Adjacency matrix


  freq = zeros();

  for i = 1:maxiter
    p = randperm(N);                  % Get permutation
    rc = r(p,p);                      % Permute matrix
    obj = corrobj(rc,w);              % Objective function for permuted matrix

    % Index into function-1 bin
    f1 = 
    distrib(d) = distrib(d)+1;        % Increment bin

  end;
  return;
