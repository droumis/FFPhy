% PARTIAL: Part of a Wright-style factor analysis.
%          Extracts 1 or 2 factors (sets of loadings) from a correlation
%          matrix, using the method of Wright (1954).
%
%     Syntax: rc = partial(c,factor1,{factor2})
%
%          c =       [n x n] correlation matrix
%          factor1 = 1st vector of loadings of length n
%          factor2 = optional 2nd vector of loadings of length n
%          -----------------------------------------------------
%          rc =      [n x n] residual correlation matrix
%

% RE Strauss, 12/23/95

function [rc,fcorr] = partial(c,factor1,factor2)
  [n,nvar] = size(c);
  if (n ~= nvar)
    error('  Error: c must be a correlation matrix');
  end;
  fcorr = [];

  if (nargin < 3)                     % Partial one factor
    for i = 1:nvar
      for j = 1:i
        rc(i,j) = c(i,j) - (factor1(i)*factor1(j));
        rc(j,i) = rc(i,j);
      end;
    end;

  else                                  % Partial two factors
    s11 = 0;
    s12 = 0;
    s22 = 0;

    % Correlation between factors
    for i = 1:nvar
      for j = 1:nvar
        cell = c(i,j);
        s11 = s11 + (cell*factor1(i)*factor1(j));
        s12 = s12 + (cell*factor1(i)*factor2(j));
        s22 = s22 + (cell*factor2(i)*factor2(j));
      end;
    end;

    if (s11*s22 > 0)
      c12 = s12/sqrt(s11*s22);
    else
      c12 = 0;
    end;
    fcorr = c12;

    % Partial the factors
    for i = 1:nvar
      for j = 1:i
        rc(i,j) = c(i,j) - (factor1(i)*c12*factor2(j)) ...
                         - (factor1(j)*c12*factor2(i));
        rc(j,i) = rc(i,j);
      end;
    end;
  end;

  return;

