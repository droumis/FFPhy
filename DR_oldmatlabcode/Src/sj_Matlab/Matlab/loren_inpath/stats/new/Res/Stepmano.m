% Stepmano: MANOVA probabilities from the results of a stepwise discriminant analysis.
%
%     Usage: pr = stepmano(X,grps,incl)
%
%         X =    [n x p] data matrix (obs x vars).
%         grps = [n x 1] vector of group identifiers, or [n x 2] matrix of group 
%                  (col 1) and subgroup (col 2) identifiers. 
%         incl = list of indices of variables in order of inclusion.
%         -----------------------------------------------------------------------------
%         pr =   cumulative MANOVA probabilities corresponding to the elements of 'incl'.
%

% RE Strauss, 6/30/98
%   11/29/99 - changed calling sequence.
%    6/21/00 - added check for missing data.

function pr = stepmano(X,grps,incl)
  if (misscheck(X,grps,incl))
    error('  STEPMANO: one or more input matrices contain missing data.');
  end;

  len_incl = length(incl);
  pr = zeros(1,len_incl);

  [F,p] = anova(X(:,incl(1)),grps);
  pr(1) = p;
  
  for i = 2:len_incl
    [lambda,F,p] = manova(X(:,incl(2:i)),grps);
    pr(i) = p;
  end;

  return;
