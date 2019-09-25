% FACTORW: Final optimization of complete model, given the primary-factor
%          loadings, secondary-factor structure and loadings, orthogonality
%          criteria, and covariance matrix.  Returns the final residual
%          covariance matrix.
%          Error handling is done by function WRIGHT.
%
%     Syntax: [fact,rc] = factorw(c,fact,secnd,orthog)
%
%         c =       original [n x n] covariance or correlation matrix.
%         fact =    [n x s+1] matrix of factor loadings (primary + secondary).
%         secnd =     [n x s] matrix of boolean values (T/F = 1/0)
%                     indicating by 1's the submatrices of variables
%                     to be included in the secondary factors (s of them).
%         orthog =  flag indicating whether secondary factors are to be:
%                     0 - oblique to primary factor and one another;
%                     1 - orthogonal to primary factor but oblique to other
%                           secondary factors;
%                     2 - orthogonal to primary factor and all other secondary
%                           factors;
%                     [default = 2].
%         ----------------------------------------------------------------
%         fact =    [n x s+1] matrix of factor loadings (primary + secondary).
%         rc =      [n x n] matrix of residual covariances.
%

function [fact,rc] = factorw(c,fact,secnd,orthog)
  [nvar,nsec] = size(secnd);

  % Create parameter vector of primary + secondary loadings

  f = fact(:,1);                % Primary loadings
  for s = 1:nsec                % Append secondary loadings
    i = find(secnd(:,s));
    f = [f; fact(i,s+1)];
  end;

  % Optimize entire model

  f = fmins('factorwf',f,[],[],c,secnd,orthog);

  % Reconstruct factor matrix

  fact(:,1) = f(1:nvar);
  low = nvar+1;
  for s = 1:nsec
    i = find(secnd(:,s));
    high = low + length(i) -1;
    fact(i,s+1) = f(low:high);
    low = high+1;
  end;

  % Final residual covariance matrix

  pf = fact(:,1);
  rc = c - pf*pf';

  r2 = vectcorr(pf,fact(:,2:nsec+1)).^2; % Independent var of secondary factors
  for s = 1:nsec
    f = fact(:,s+1);
    rc = rc - (f*f')*(1-r2(s));        % Residual covariances
  end;

  return;

