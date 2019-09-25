% MEDIANS:  median value.  For column vectors, medians(x) returns the median 
%           value.  For matrices or row vectors, medians(x) is a row vector 
%           containing  the median value of each column.  The basic difference 
%           from the Matlab function median() is for a row vector, where 
%           medians() returns the row vector instead of the median value of the 
%           elements of the row.  Allows for missing data.
%
%           If an optional grouping vector is supplied, returns a vector of 
%           medians for each group.
%
%     Usage: [M,grpids] = medians(X,{grps})
%
%         X =       [n x p] data matrix.
%         grps =    optional [n x 1] grouping vector for k groups.
%         ----------------------------------------------------------------
%         M =       [k x p] matrix of group medians.
%         grpids =  corresponding [k x 1] vector of group identifiers, for 
%                     multiple groups.
%

% RE Strauss, 12/21/99

function [M,grpids] = medians(X,grps)
  if (nargin < 2) grps = []; end;

  [n,p] = size(X);

  if (isempty(grps))
    grps = ones(n,1);
    ug = 1;
    k = 1;
  else
    ug = uniquef(grps);  
    k = length(ug);
  end;

  grpids = zeros(k,1);
  M = zeros(k,p);

  for ik = 1:k
    ir = find(grps==ug(ik));
    for c = 1:p
      x = X(ir,c);
      ic = find(isfinite(x));
      if (isempty(ic))
        M(ik,c) = NaN;
      else
        M(ik,c) = median(x(ic));
      end;
    end;
  end;

  return;
