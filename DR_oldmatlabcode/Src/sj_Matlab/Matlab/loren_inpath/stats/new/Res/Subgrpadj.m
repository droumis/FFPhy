% SUBGRPADJ:  Adjust for subgroups within groups by setting subgroup means equal 
%             to the grand mean, within groups.  Adjusts each column separately.
%
%     Usage: Y = subgrpadj(X,grps,subgrps)
%
%         X =       [n x p] data matrix.
%         grps =    [n x 1] vector of group identifiers.
%         subgrps = [n x 1] vector of subgroup identifiers.
%         -------------------------------------------------
%         Y =       [n x p] adjusted data matrix.
%

% RE Strauss, 10/19/00

function Y = subgrpadj(X,grps,subgrps)
  if (nargin < 3)
    error('  SUBGRPADJ: all input matrices must be provided.');
  end;
  if (~isvector(grps) | ~isvector(subgrps))
    error('  SUBGRPADJ: group- and subgroup-identifiers must be vectors');
  end;
  
  if (isvector(X))
    X = X(:);
  end;
  [n,p] = size(X);

  if (length(grps)~=n | length(subgrps)~=n)
    error('  SUBGRPADJ: all input matrices not of same length.');
  end;

  Y = X;                                  % Allocate output matrix
  ug = uniquef(grps);                     % Group identifiers

  for ig = 1:length(ug)                   % For each group,
    i = find(grps == ug(ig));               % Isolate observations
    x = X(i,:);
    sg = subgrps(i);                        % Identify subgroups
    xm = means(x);                          % Grand means
    
    usg = uniquef(sg);
    for isg = 1:length(usg)                 % For each subgroup,
      j = find(sg == usg(isg));               % Identify observations
      dm = means(x(j,:)) - xm;                % Subgroup delta-means
      x(j,:) = x(j,:) - dm;                   % Adjust subgroup means
    end;
    Y(i,:) = x;
  end;

  return;
