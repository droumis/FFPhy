% SUMS:   Column sums.  For column vectors, sums(x) returns the 
%         sum value.  For matrices or row vectors, sums(x) is a row vector 
%         containing the sum of each column.  The basic difference from 
%         the Matlab function sum() is for a row vector, where sums() returns 
%         the row vector instead of the sum value of the elements of the row.   
%         Allows for missing data.
%
%         If an optional grouping vector is supplied, returns a vector of sums
%         for each group in collating sequence.
%
%     Usage: [S,grpids] = sums(X,{grps})
%
%         X =       [n x p] data matrix.
%         grps =    optional [n x 1] grouping vector for k groups.
%         ---------------------------------------------------------------------
%         S =       [k x p] matrix of group sums.
%         grpids =  corresponding [k x 1] vector of group identifiers, for multiple
%                     groups.
%

% RE Strauss, 10/20/00 - modified from means().

function [S,grpids] = sums(X,grps)
  if (nargin < 2) grps = []; end;

  [n,p] = size(X);

  if (isempty(grps))
    grps = ones(n,1);
    ug = 1;
    k = 1;
  else
    [ug,fg] = uniquef(grps,1);  
    k = length(ug);
  end;

  grpids = ug;
  S = zeros(k,p);

  for ik = 1:k
    ir = find(grps==ug(ik));
    for c = 1:p
      x = X(ir,c);
      ic = find(isfinite(x));
      if (isempty(ic))
        S(ik,c) = NaN;
      else
        S(ik,c) = sum(x(ic));
      end;
    end;
  end;

  return;
