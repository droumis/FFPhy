% RANDPERMG: Randomly permute the rows of a matrix, optionally by group.
%
%     Usage: Y = randpermg(X,{grps})
%
%         X =     [n x p] data matrix.
%         grps =  optional [n x 1] group-identification vector.
%         -------------------------------------------------
%         Y =     permuted matrix.
%

% RE Strauss, 12/3/00

function Y = randpermg(X,grps)
  if (nargin < 2) grps = []; end;

  [n,p] = size(X);
  if (isempty(grps))
    grps = ones(n,1);
  end;

  if (~isvector(grps))
    error('  RANDPERMG: group-id matrix must be a vector.');
  end;
  if (length(grps)~=n)
    error('  RANDPERMG: input matrices not compatible in size.');
  end;

  ugrps = uniquef(grps);
  ngrps = length(ugrps);

  for g = 1:ngrps                     % Cycle thru groups
    i = find(grps == ugrps(g));         % Isolate current group
    x = X(i,:);
    p = randperm(length(i));
    Y(i,:) = x(p,:);                    % Permute into output matrix
  end;

  return;

    